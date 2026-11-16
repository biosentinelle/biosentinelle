# Web page of biosentinelle: Shiny form (web/app/app.R) + submission daemon to Galaxy (web/daemon/bs_bioblend.py).
# Build from the root of the repository:
#   docker build -f web/Dockerfile -t biosentinelle_web .
# Run:
#   docker run -p 3838:3838 -e GALAXY_API_KEY=... -e BS_PUBLIC_URL=https://... -v bs_work:/srv/bs_work biosentinelle_web
FROM rocker/shiny:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --break-system-packages bioblend==1.9.0 \
    && rm -rf /srv/shiny-server/*

# only the files used by the web page and the daemon (the pipeline itself runs in Galaxy)
COPY nextflow.config /opt/biosentinelle/nextflow.config
COPY galaxy/fill_config.py galaxy/locked_params.txt /opt/biosentinelle/galaxy/
COPY web/daemon/bs_bioblend.py /opt/biosentinelle/web/daemon/
COPY web/app/app.R /srv/shiny-server/app.R
COPY web/start.sh /usr/bin/start.sh

RUN chmod +x /usr/bin/start.sh \
    && mkdir -p /srv/bs_work && chown -R shiny:shiny /srv/bs_work /srv/shiny-server

# BS_MAX_UPLOAD_MB: total size of the sample files of a submission (sequencing runs are large: the run accession
# option of the page avoids the upload)
ENV BS_WORK_DIR=/srv/bs_work \
    BS_REPO_DIR=/opt/biosentinelle \
    BS_PYTHON=python3 \
    BS_MAX_UPLOAD_MB=1000 \
    BS_RETENTION_DAYS=15

VOLUME /srv/bs_work
EXPOSE 3838
CMD ["/usr/bin/start.sh"]
