#!/usr/bin/env python3
"""
Progress of the nextflow processes of the biosentinelle jobs run by the local Galaxy (serve.sh, test.sh),
in a terminal, like the console of `nextflow run main.nf`.

The Galaxy tool runs nextflow with -ansi-log false, so that its standard output stays readable in the Galaxy page
(one line per started task, shown live in the dataset details thanks to live_tool_output_reporting in settings.sh).
This script rebuilds the progress table of the nextflow console from the .nextflow.log file of the job, every 2 seconds:

    executor >  local (25)
    [3f/2a1b3c] PrepareReads (1)       | 1 of 1 ✔
    [a2/b3c4d5] Bowtie2_ref (2)        | 1 of 2
    [-        ] StrainScore            -

followed by the last warnings and errors of nextflow. Differences with the nextflow console:
    - "of N" counts the tasks already submitted (nextflow also counts the tasks waiting for free cores)
    - ✔ means that all the submitted tasks of the process are done (others can come until the end of the run)
At the end of the job, Galaxy deletes the job folder: the last table stays on the screen and the script waits for the next job.

Usage, in a second terminal, from any folder (Ctrl+C to stop):
    python3 galaxy/local/watch.py            # last started job, then the next ones
    python3 galaxy/local/watch.py <folder>   # only the job of this Galaxy job folder (or of its working subfolder)
    python3 galaxy/local/watch.py --once     # prints the table once
Settings: BS_WORK_DIR and BS_TMPDIR, as for serve.sh and test.sh (see settings.sh).
The Galaxy job folders are BS_TMPDIR/tmp*/job_working_directory/*/<job id>, with the nextflow launch folder in working/.
"""

import argparse
import glob
import os
import re
import shutil
import sys
import time

REFRESH = 2  # seconds between two refreshes of the screen
RECENT = 600  # seconds: jobs whose .nextflow.log did not change for longer are ignored (nextflow writes it at least every 5 minutes)
MAX_MESSAGES = 3  # number of warnings and errors shown under the table

# line of .nextflow.log, e.g., "Nov-15 19:11:01.013 [Task submitter] INFO  nextflow.Session - [47/a548e9] Submitted process > A (2)"
LOG_LINE = re.compile(r"^[A-Z][a-z]{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3} \[[^\]]*\] (?P<level>[A-Z]+) +\S+ - (?P<msg>.*)$")
# messages of the log used to rebuild the progress
VERSION = re.compile(r"^N E X T F L O W\s+~\s+version (?P<version>\S+)")
LAUNCH = re.compile(r"^Launching `")
EXECUTOR = re.compile(r"^\[warm up\] executor > (?P<executor>\S+)")
START = re.compile(r"^Starting process > (?P<process>\S+)")
TASK = re.compile(r"^\[(?P<hash>[0-9a-f]{2}/[0-9a-f]{6})\] (?P<kind>Submitted|Re-submitted|Cached) process > (?P<task>.+)$")
COMPLETED = re.compile(r"^Task completed > \w*TaskHandler\[(?:.*?; )?id: \d+; name: (?P<task>.+?); status: \w+; exit: (?P<exit>[^;]*);")
RETRIED = re.compile(r"Process `(?P<task>.+?)` terminated .*-- Execution is retried")
FAILED = re.compile(r"^Error executing process > '(?P<task>.+)'$")
STATS = re.compile(r"^Workflow completed > WorkflowStats\[(?P<stats>.*)\]")
GOODBYE = re.compile(r"Execution complete -- Goodbye")
# sections of the nextflow errors shown under the table (number of lines, from the end or not); full error: Galaxy page or .nextflow.log
ERROR_SECTIONS = {"Caused by:": (2, False), "Command error:": (5, True), "Work dir:": (1, False)}
STYLES = {"dim": "\033[2m", "red": "\033[31m", "yellow": "\033[33m", "green": "\033[32m"}


class Run:
    """Progress of a nextflow run, rebuilt from its .nextflow.log file, read as it grows."""

    def __init__(self, log):
        self.log = log  # path of .nextflow.log
        self.offset = 0  # bytes of the log already read
        self.rest = b""  # incomplete last line, completed by the next read
        self.version = self.launch = self.stats = None
        self.executor = "local"
        self.finished = False  # end of the nextflow run
        self.processes = {}  # process name: counts of its tasks, in the order of the nextflow console
        self.exits = {}  # task name: exit status of its last execution
        self.messages = []  # last warnings and errors: [level, lines]
        self.message = None  # message receiving the next lines of the log without date (multi-line message)

    def counts(self, task):
        """Counts of the process of a task (task name: process name, then its tag between brackets)."""
        process = task.split(" ", 1)[0]
        return self.processes.setdefault(
            process, {"hash": None, "task": process, "submitted": 0, "cached": 0, "succeeded": 0, "failed": 0, "retries": 0}
        )

    def update(self):
        """Reads the new lines of the log. Returns False if the log does not exist anymore (job folder deleted by Galaxy)."""
        try:
            with open(self.log, "rb") as handle:
                if os.fstat(handle.fileno()).st_size < self.offset:  # log written again from the start
                    self.__init__(self.log)
                handle.seek(self.offset)
                data = handle.read()
        except OSError:
            return False
        self.offset += len(data)
        lines = (self.rest + data).split(b"\n")
        self.rest = lines.pop()
        for line in lines:
            self.parse(line.decode("utf-8", "replace").rstrip("\r"))
        return True

    def parse(self, line):
        """Updates the progress with a line of the log."""
        match = LOG_LINE.match(line)
        if match is None:  # next line of a multi-line message
            if self.message is not None and len(self.message[1]) < 500:
                self.message[1].append(line)
            return
        level, msg = match["level"], match["msg"]
        self.message = None
        if m := TASK.match(msg):
            counts = self.counts(m["task"])
            counts["hash"], counts["task"] = m["hash"], m["task"]
            counts["cached" if m["kind"] == "Cached" else "submitted"] += 1
        elif m := COMPLETED.match(msg):
            self.exits[m["task"]] = m["exit"]
            self.counts(m["task"])["succeeded" if m["exit"] == "0" else "failed"] += 1
        elif m := RETRIED.search(msg):  # the failed task is submitted again
            counts = self.counts(m["task"])
            counts["failed"] = max(counts["failed"] - 1, 0)
            counts["retries"] += 1
        elif (m := FAILED.match(msg)) and self.exits.get(m["task"]) == "0":  # error after a 0 exit status, e.g., missing output file
            counts = self.counts(m["task"])
            counts["succeeded"] -= 1
            counts["failed"] += 1
            self.exits[m["task"]] = "error"
        elif m := START.match(msg):
            self.counts(m["process"])
        elif m := VERSION.match(msg):
            self.version = m["version"]
        elif LAUNCH.match(msg):
            self.launch = msg
        elif m := EXECUTOR.match(msg):
            self.executor = m["executor"]
        elif m := STATS.match(msg):
            self.stats = dict(re.findall(r"(\w+)=([^;]*)", m["stats"]))
        elif GOODBYE.search(msg):
            self.finished = True
        if level in ("WARN", "ERROR"):
            self.message = [level, [msg]]
            self.messages = (self.messages + [self.message])[-MAX_MESSAGES:]


def summary(level, lines):
    """Short version of a nextflow warning or error: its first line, then the main sections of the process errors."""
    out = [("ERROR ~ " if level == "ERROR" else "WARN: ") + lines[0]]
    sections, name = {}, None
    for line in lines[1:]:
        if line and not line[0].isspace() and line.endswith(":"):  # section title, e.g., "Caused by:"
            name = line
            sections[name] = []
        elif line.strip():
            sections.setdefault(name, []).append(line.strip())
    if any(key in sections for key in ERROR_SECTIONS):
        for key, (number, from_end) in ERROR_SECTIONS.items():
            kept = sections.get(key, [])
            kept = kept[-number:] if from_end else kept[:number]
            if kept:
                out += [f"  {key} {kept[0]}"] + [f"    {line}" for line in kept[1:]]
    else:
        out += [f"  {line}" for line in sections.get(None, [])[:3]]
    return out


def process_lines(run, columns):
    """Lines of the processes, as in the nextflow console, with their state: waiting, running, done or failed."""
    width = min(max((len(counts["task"]) for counts in run.processes.values()), default=0), max(columns - 40, 15))
    lines = []
    for counts in run.processes.values():
        label = counts["task"] if len(counts["task"]) <= width else counts["task"][: width - 1] + "…"
        if counts["hash"] is None:
            lines.append((f"[-        ] {label:<{width}} -", "waiting"))
            continue
        total = counts["submitted"] + counts["cached"] - counts["retries"]
        done = counts["succeeded"] + counts["cached"] + counts["failed"]
        text = f"[{counts['hash']}] {label:<{width}} | {done} of {total}"
        text += "".join(f", {key}: {counts[key]}" for key in ("cached", "failed", "retries") if counts[key])
        if counts["failed"]:
            lines.append((text + " ✘", "failed"))
        elif done >= total:
            lines.append((text + " ✔", "done"))
        else:
            lines.append((text, "running"))
    return lines


def fit(lines, room):
    """Process lines reduced to room lines: the processes not started, then those done, are summarized on one line."""
    if len(lines) > room:
        waiting = [line for line in lines if line[1] == "waiting"]
        lines = [line for line in lines if line[1] != "waiting"]
        if waiting:
            lines.append((f"[-        ] {len(waiting)} other process{'es' if len(waiting) > 1 else ''} not started", "waiting"))
    if len(lines) > room:
        done = [line for line in lines if line[1] == "done"]
        lines = [line for line in lines if line[1] != "done"]
        if done:
            lines.insert(0, (f"[         ] {len(done)} other process{'es' if len(done) > 1 else ''} done ✔", "done"))
    if len(lines) > room:
        lines = lines[: max(room - 1, 0)] + [("... (enlarge the terminal to see the other processes)", "waiting")]
    return lines


def screen(run, title, status, columns, rows):
    """Lines of the screen, as (text, style)."""
    head, body, tail = [(title, "dim"), ("", None)], [], []
    if run is not None:
        if run.version:
            head += [(f" N E X T F L O W   ~  version {run.version}", None), ("", None)]
        if run.launch:
            head += [(run.launch, None), ("", None)]
        if run.processes:
            head.append((f"executor >  {run.executor} ({sum(counts['submitted'] for counts in run.processes.values())})", None))
        body = process_lines(run, columns)
        for level, lines in run.messages:
            tail += [("", None)] + [(line, "red" if level == "ERROR" else "yellow") for line in summary(level, lines)]
        if run.stats:  # end of the run: counts of nextflow, whose failed tasks include the failed attempts of the retried tasks
            stats = {key: int(value) for key, value in run.stats.items() if value.isdigit()}
            error = any(level == "ERROR" for level, _ in run.messages) or any(c["failed"] for c in run.processes.values())
            text = f"Nextflow run {'ended with an error' if error else 'completed'}: {stats.get('succeededCount', 0)} tasks succeeded"
            numbers = {
                "cached": stats.get("cachedCount", 0),
                "failed": max(stats.get("failedCount", 0) - stats.get("retriesCount", 0), 0),
                "retried after a failure": stats.get("retriesCount", 0),
                "ignored": stats.get("ignoredCount", 0),
            }
            text += "".join(f", {number} {name}" for name, number in numbers.items() if number)
            tail += [("", None), (text, "red" if error else "green")]
    tail += [("", None), (status, "dim")]
    return head + fit(body, rows - 1 - len(head) - len(tail)) + tail


def draw(lines, columns, tty):
    """Prints the screen: drawn again in place in a terminal, printed after the previous one otherwise."""
    if not tty:
        sys.stdout.write("\n".join(text for text, _ in lines) + "\n\n")
        sys.stdout.flush()
        return
    out = []
    for text, style in lines:
        text = text[:columns]
        if style in STYLES:
            text = STYLES[style] + text + "\033[0m"
        else:
            text = text.replace(" ✔", " \033[32m✔\033[0m").replace(" ✘", " \033[31m✘\033[0m")
        out.append(text + "\033[K")  # \033[K: end of the previous screen line erased
    sys.stdout.write("\033[H" + "\n".join(out) + "\n\033[J")  # \033[H: top of the screen, \033[J: rest of the screen erased
    sys.stdout.flush()


def resolve(folder):
    """.nextflow.log of a Galaxy job folder, of its working subfolder, or the file itself (None if absent)."""
    for path in (os.path.join(folder, "working", ".nextflow.log"), os.path.join(folder, ".nextflow.log"), folder):
        if os.path.basename(path) == ".nextflow.log" and os.path.isfile(path):
            return os.path.abspath(path)
    return None


def newest(pattern, excluded):
    """.nextflow.log written last among the recent ones not followed yet (None if none)."""
    now, found = time.time(), []
    for log in glob.glob(pattern):
        try:
            mtime = os.path.getmtime(log)
        except OSError:  # job folder deleted in between
            continue
        if log not in excluded and now - mtime < RECENT:
            found.append((mtime, log))
    return max(found)[1] if found else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", nargs="?", help="Galaxy job folder to follow (default: the last started job, then the next ones)")
    parser.add_argument("--once", action="store_true", help="print the progress once, then stop")
    args = parser.parse_args()

    work_dir = os.environ.get("BS_WORK_DIR") or os.path.join(os.path.expanduser("~"), ".planemo", "biosentinelle")
    tmp_dir = os.environ.get("BS_TMPDIR") or os.path.join(work_dir, "tmp")
    pattern = os.path.join(tmp_dir, "tmp*", "job_working_directory", "*", "*", "working", ".nextflow.log")
    pinned = None
    if args.folder:
        pinned = resolve(args.folder)
        if pinned is None:
            sys.exit(f"ERROR: no .nextflow.log in {args.folder} or in its working subfolder")

    tty = sys.stdout.isatty() and not args.once
    run, gone, followed = None, False, set()  # gone: job folder deleted by Galaxy; followed: logs of the jobs already ended
    if tty:
        sys.stdout.write("\033[?25l\033[2J")  # cursor hidden, screen erased
    try:
        while True:
            if pinned:
                run = run or Run(pinned)
            elif run is None or gone or run.finished:  # a new job replaces the one displayed
                log = newest(pattern, followed | ({run.log} if run else set()))
                if log:
                    run, gone = Run(log), False
            if run is not None and not gone:
                gone = not run.update()
                if gone or run.finished:
                    followed.add(run.log)
            columns, rows = shutil.get_terminal_size((120, 50))
            if run is None:
                title = "biosentinelle jobs of the local Galaxy: progress of the nextflow processes"
                status = f"No running job in {os.path.dirname(os.path.dirname(pattern))}. Waiting for a job (Ctrl+C to stop)"
            else:
                launch_dir = os.path.dirname(run.log)
                if os.path.basename(launch_dir) == "working":  # Galaxy job folder: <job id>/working/.nextflow.log
                    title = f"Galaxy job {os.path.basename(os.path.dirname(launch_dir))} ({os.path.dirname(launch_dir)})"
                else:
                    title = f"Nextflow run of {launch_dir}"
                title += f"  {time.strftime('%H:%M:%S')}"
                if gone:
                    status = "Galaxy job ended (job folder deleted by Galaxy): see the result in the Galaxy history."
                    status += "" if pinned else " Waiting for the next job (Ctrl+C to stop)"
                elif run.finished:
                    status = "Nextflow run finished (the Galaxy job then copies the results)"
                else:
                    status = f"Refreshed every {REFRESH} seconds (Ctrl+C to stop)"
            draw(screen(run, title, status, columns, rows), columns, tty)
            if args.once or (pinned and (gone or run.finished)):
                break
            time.sleep(REFRESH)
    except KeyboardInterrupt:
        pass
    finally:
        if tty:
            sys.stdout.write("\033[?25h")  # cursor shown again
            sys.stdout.flush()


if __name__ == "__main__":
    main()
