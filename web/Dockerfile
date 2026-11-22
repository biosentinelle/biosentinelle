# Web page of lab_strain_watch: Shiny form (web/app/app.R) + submission daemon to Galaxy (web/daemon/lsw_bioblend.py).
# Build from the root of the repository:
#   docker build -f web/Dockerfile -t lab_strain_watch_web .
# Run:
#   docker run -p 3838:3838 -e GALAXY_API_KEY=... -e LSW_PUBLIC_URL=https://... -v lsw_work:/srv/lsw_work lab_strain_watch_web
FROM rocker/shiny:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --break-system-packages bioblend==1.9.0 \
    && rm -rf /srv/shiny-server/*

# only the files used by the web page and the daemon (the pipeline itself runs in Galaxy)
COPY nextflow.config /opt/lab_strain_watch/nextflow.config
COPY galaxy/fill_config.py galaxy/locked_params.txt /opt/lab_strain_watch/galaxy/
COPY web/daemon/lsw_bioblend.py /opt/lab_strain_watch/web/daemon/
COPY web/app/app.R /srv/shiny-server/app.R
COPY web/start.sh /usr/bin/start.sh

RUN chmod +x /usr/bin/start.sh \
    && mkdir -p /srv/lsw_work && chown -R shiny:shiny /srv/lsw_work /srv/shiny-server

# LSW_MAX_UPLOAD_MB: total size of the sample files of a submission (sequencing runs are large: the run accession
# option of the page avoids the upload)
ENV LSW_WORK_DIR=/srv/lsw_work \
    LSW_REPO_DIR=/opt/lab_strain_watch \
    LSW_PYTHON=python3 \
    LSW_MAX_UPLOAD_MB=1000 \
    LSW_RETENTION_DAYS=15

VOLUME /srv/lsw_work
EXPOSE 3838
CMD ["/usr/bin/start.sh"]
