#!/usr/bin/env python3
"""
Submission daemon of the lab_strain_watch web page (same principle as shaman_bioblend.py of SHAMAN).

The web page (web/app/app.R) writes, for each submission with key KEY:
    <work_dir>/jobs/KEY/sample/*          sequence files of the sample (absent with a run accession)
    <work_dir>/jobs/KEY/params.json       values of the form
    <work_dir>/todo/KEY.json              task: {"key", "mail", "accession", "sample_name", "submitted"}

This daemon polls todo/, moves each task to doing/, and in a thread:
    1) creates a Galaxy history and uploads the files, under their names (the mates are paired by their _1/_2 or
       _R1/_R2 suffix),
    2) runs the lab_strain_watch Galaxy tool (galaxy/lab_strain_watch.xml): the files (or the run accession) are its sample
       input and params.json its parameter file (whose values replace those of the tool form, whose defaults are those
       of nextflow.config),
    3) waits for the job, downloads the results into jobs/KEY/ (results.zip, the result folder, and
       strain_identity_report.html and report.html, extracted from it) and moves the task to done/ (or error/),
    4) sends an e-mail if an address was given.
The task file always holds the current "state" and "message", read by the web page for the "Check my job" tab.
Jobs older than LSW_RETENTION_DAYS are removed, locally and in Galaxy.

Configuration (environment variables, never in the code):
    GALAXY_API_KEY      required, key of the Galaxy service account
    GALAXY_URL          default https://galaxy.pasteur.fr
    LSW_TOOL_ID         default lab_strain_watch (full Tool Shed id if installed from a Tool Shed)
    LSW_MAX_JOBS        max number of jobs running at the same time in Galaxy, default 5
    LSW_RETENTION_DAYS  default 15
    LSW_PUBLIC_URL      URL of the web page, written in the e-mails
    LSW_SMTP_HOST       default smtp.pasteur.fr ("" for no e-mail)
    LSW_MAIL_FROM       default lab_strain_watch@pasteur.fr

Usage:
    lsw_bioblend.py -w <work_dir> [--once]
"""

import argparse
import datetime
import json
import logging
import os
import re
import shutil
import smtplib
import sys
import threading
import time
import zipfile
from email.message import EmailMessage

import requests
from bioblend import ConnectionError as GalaxyConnectionError
from bioblend.galaxy import GalaxyInstance

KEY_RE = re.compile(r"^[A-Za-z0-9_-]{16,64}$")
ACCESSION_RE = re.compile(r"^(SRR|ERR|DRR)[0-9]+$")
POLL_SECONDS = 5
JOB_POLL_SECONDS = 30
RETRY_SECONDS = 300  # when Galaxy is unreachable (e.g., maintenance)
UPLOAD_MAX_WAIT = 6 * 3600  # seconds: sequencing runs can be large
REPORTS = ("strain_identity_report.html", "report.html")  # extracted from the result zip for the web page

log = logging.getLogger("lsw_bioblend")


class Config:
    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        self.galaxy_url = os.environ.get("GALAXY_URL", "https://galaxy.pasteur.fr")
        self.galaxy_key = os.environ.get("GALAXY_API_KEY", "")
        self.tool_id = os.environ.get("LSW_TOOL_ID", "lab_strain_watch")
        self.max_jobs = int(os.environ.get("LSW_MAX_JOBS", "5"))
        self.retention_days = int(os.environ.get("LSW_RETENTION_DAYS", "15"))
        self.public_url = os.environ.get("LSW_PUBLIC_URL", "")
        self.smtp_host = os.environ.get("LSW_SMTP_HOST", "smtp.pasteur.fr")
        self.mail_from = os.environ.get("LSW_MAIL_FROM", "lab_strain_watch@pasteur.fr")
        if not self.galaxy_key:
            sys.exit("lsw_bioblend.py: ERROR: the GALAXY_API_KEY environment variable is not set")

    def path(self, *parts):
        return os.path.join(self.work_dir, *parts)


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def read_task(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_task(path, task):
    # write then rename, so that the web page never reads a half written file
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(task, f, indent=1)
    os.replace(tmp, path)


class Job(threading.Thread):
    """One submission, from doing/KEY.json to done/ or error/."""

    def __init__(self, cfg, key, slots):
        super().__init__(name="job_" + key, daemon=True)
        self.cfg = cfg
        self.key = key
        self.slots = slots
        self.task_file = cfg.path("doing", key + ".json")
        self.job_dir = cfg.path("jobs", key)
        self.task = read_task(self.task_file)

    def update(self, state, message):
        self.task.update({"state": state, "message": message, "updated": now()})
        write_task(self.task_file, self.task)
        log.info("%s: %s: %s", self.key, state, message)

    def finish(self, folder, state, message):
        self.update(state, message)
        dest = self.cfg.path(folder, self.key + ".json")
        os.replace(self.task_file, dest)
        self.task_file = dest
        self.send_mail(state, message)

    def run(self):
        with self.slots:
            while True:
                try:
                    self.process()
                    return
                except (GalaxyConnectionError, requests.exceptions.RequestException) as e:
                    status = getattr(e, "status_code", None)
                    if status is not None and status < 500:
                        log.exception("%s: Galaxy refused the request", self.key)
                        self.finish("error", "error", "Galaxy refused the job (HTTP {0}). Please contact the administrators.".format(status))
                        return
                    # Galaxy down or under maintenance: keep the task in doing/ and try again later
                    log.warning("%s: Galaxy connection error, retry in %s s: %s", self.key, RETRY_SECONDS, e)
                    self.update("waiting", "The Galaxy server is not reachable for the moment. The job will be submitted again automatically.")
                    time.sleep(RETRY_SECONDS)
                except Exception as e:
                    log.exception("%s: failure", self.key)
                    self.finish("error", "error", "Internal error of the submission system: {0}".format(e))
                    return

    def process(self):
        gi = GalaxyInstance(url=self.cfg.galaxy_url, key=self.cfg.galaxy_key)
        if not self.task.get("job_id"):
            self.submit(gi)
        self.wait(gi)

    def upload(self, gi, hid, path, file_type):
        out = gi.tools.upload_file(path, hid, file_type=file_type, file_name=os.path.basename(path))
        return {"src": "hda", "id": out["outputs"][0]["id"]}

    def submit(self, gi):
        if not self.task.get("history_id"):
            history = gi.histories.create_history(name="lab_strain_watch_" + self.key)
            self.task["history_id"] = history["id"]
            self.update("uploading", "Uploading the files to Galaxy.")
        hid = self.task["history_id"]
        accession = self.task.get("accession", "")
        sample_dir = os.path.join(self.job_dir, "sample")
        # inputs of the tool form (sections and conditionals of galaxy/lab_strain_watch.xml); the other fields keep
        # their default (that of nextflow.config) and the parameter file replaces them with the values of the web page
        if accession:
            if not ACCESSION_RE.match(accession):
                raise ValueError("invalid run accession {0!r}".format(accession))
            inputs = {"input": {"input_type": "accession", "accession": accession}}
            datasets = []
        else:
            # zip archives as they are, the other files with the datatype detected by Galaxy (e.g., fastqsanger.gz)
            names = sorted(os.listdir(sample_dir))
            if not names:
                raise ValueError("no sample file in {0}".format(sample_dir))
            datasets = [self.upload(gi, hid, os.path.join(sample_dir, n), "zip" if n.lower().endswith(".zip") else "auto") for n in names]
            inputs = {"input": {"input_type": "files", "seq_files": {"values": datasets}, "sample_name": self.task.get("sample_name", "")}}
        inputs["output_options"] = {"zip": True, "tables": False}
        params_file = os.path.join(self.job_dir, "params.json")
        if os.path.isfile(params_file):
            params = self.upload(gi, hid, params_file, "json")
            datasets.append(params)
            inputs["advanced"] = {"params": params}
        for value in datasets:
            gi.datasets.wait_for_dataset(value["id"], maxwait=UPLOAD_MAX_WAIT)
        run = gi.tools.run_tool(hid, self.cfg.tool_id, inputs, input_format="21.01")
        self.task["job_id"] = run["jobs"][0]["id"]
        self.task["outputs"] = {o["output_name"]: o["id"] for o in run["outputs"]}
        self.update("running", "The analysis is running on the Institut Pasteur cluster.")

    def wait(self, gi):
        job_id = self.task["job_id"]
        while True:
            job = gi.jobs.show_job(job_id)
            state = job["state"]
            if state == "ok":
                break
            if state in ("error", "deleted", "deleting", "failed", "stopped"):
                details = gi.jobs.show_job(job_id, full_details=True)
                with open(os.path.join(self.job_dir, "error.txt"), "w", encoding="utf-8") as f:
                    f.write("stdout:\n{0}\n\nstderr:\n{1}\n".format(details.get("stdout", ""), details.get("stderr", "")))
                self.finish("error", "error", "The analysis failed. See the error log of the job, and check your files and parameters.")
                return
            if self.task.get("state") != "running":
                self.update("running", "The analysis is running on the Institut Pasteur cluster.")
            time.sleep(JOB_POLL_SECONDS)
        self.update("downloading", "Downloading the results.")
        zip_path = os.path.join(self.job_dir, "results.zip")
        gi.datasets.download_dataset(self.task["outputs"]["results_zip"], file_path=zip_path, use_default_filename=False)
        # the zip holds the result folder lab_strain_watch_<time>/, with the two HTML reports at its root
        with zipfile.ZipFile(zip_path) as z:
            for report in REPORTS:
                members = [n for n in z.namelist() if re.fullmatch(r"[^/]+/" + re.escape(report), n)]
                if members:
                    with z.open(members[0]) as src, open(os.path.join(self.job_dir, report), "wb") as dest:
                        shutil.copyfileobj(src, dest)
        self.finish("done", "done", "The analysis is finished.")

    def send_mail(self, state, message):
        mail = self.task.get("mail", "")
        if not mail or not self.cfg.smtp_host:
            return
        msg = EmailMessage()
        msg["Subject"] = "lab_strain_watch job {0}: {1}".format(self.key, state)
        msg["From"] = self.cfg.mail_from
        msg["To"] = mail
        text = "{0}\n\nYour job key is: {1}\n".format(message, self.key)
        if self.cfg.public_url:
            text += "Results and status: {0}?key={1}\n".format(self.cfg.public_url, self.key)
        text += "\nResults are kept {0} days.\n".format(self.cfg.retention_days)
        msg.set_content(text)
        try:
            with smtplib.SMTP(self.cfg.smtp_host, 25, timeout=60) as server:
                server.send_message(msg)
        except Exception:
            log.exception("%s: cannot send the e-mail", self.key)


def cleanup(cfg, gi):
    """Remove the jobs (local files and Galaxy history) older than the retention time."""
    limit = time.time() - cfg.retention_days * 86400
    for folder in ("done", "error"):
        for name in os.listdir(cfg.path(folder)):
            if not name.endswith(".json"):
                continue
            path = cfg.path(folder, name)
            if os.path.getmtime(path) > limit:
                continue
            key = name[:-5]
            try:
                task = read_task(path)
                if task.get("history_id"):
                    gi.histories.delete_history(task["history_id"], purge=True)
            except Exception:
                log.exception("%s: cannot delete the Galaxy history", key)
                continue
            shutil.rmtree(cfg.path("jobs", key), ignore_errors=True)
            os.remove(path)
            log.info("%s: removed (retention time)", key)


def start_task(cfg, key, slots, running):
    if key in running and running[key].is_alive():
        return
    job = Job(cfg, key, slots)
    running[key] = job
    job.start()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-w", dest="work_dir", required=True, help="top folder shared with the web page")
    parser.add_argument("--once", action="store_true", help="process the waiting tasks, wait for them and exit (tests)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(threadName)s %(message)s")
    cfg = Config(args.work_dir)
    for folder in ("todo", "doing", "done", "error", "jobs"):
        os.makedirs(cfg.path(folder), exist_ok=True)

    slots = threading.BoundedSemaphore(cfg.max_jobs)
    running = {}
    # tasks interrupted by a restart of the daemon: resume them
    for name in sorted(os.listdir(cfg.path("doing"))):
        if name.endswith(".json"):
            start_task(cfg, name[:-5], slots, running)

    last_cleanup = 0
    while True:
        for name in sorted(os.listdir(cfg.path("todo"))):
            if not name.endswith(".json"):
                continue
            key = name[:-5]
            src = cfg.path("todo", name)
            if not KEY_RE.match(key) or not os.path.isdir(cfg.path("jobs", key)):
                log.error("invalid task %s, moved to error/", name)
                os.replace(src, cfg.path("error", name))
                continue
            dest = cfg.path("doing", name)
            os.replace(src, dest)
            task = read_task(dest)
            task.update({"state": "queued", "message": "The job is waiting to be sent to Galaxy.", "updated": now()})
            write_task(dest, task)
            start_task(cfg, key, slots, running)
        if time.time() - last_cleanup > 3600:
            try:
                cleanup(cfg, GalaxyInstance(url=cfg.galaxy_url, key=cfg.galaxy_key))
            except Exception:
                log.exception("cleanup failed")
            last_cleanup = time.time()
        if args.once and not any(j.is_alive() for j in running.values()):
            return
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    main()
