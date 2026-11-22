#!/bin/sh
# Starts the submission daemon and Shiny Server (as SHAMAN's docker_inst/shiny-server.sh).
set -e

if [ -z "$GALAXY_API_KEY" ]; then
    echo "start.sh: ERROR: the GALAXY_API_KEY environment variable is not set" >&2
    exit 1
fi

# Shiny Server does not pass the container environment to the R processes: write the settings of app.R (not the Galaxy key)
env | grep '^LSW_' >> "$(R RHOME)/etc/Renviron.site"

chown -R shiny:shiny "$LSW_WORK_DIR"

# daemon as the shiny user, so that both processes can read and write the work folder. Restarted if it stops
(
    while true; do
        runuser -u shiny -- python3 "$LSW_REPO_DIR/web/daemon/lsw_bioblend.py" -w "$LSW_WORK_DIR" >> /var/log/shiny-server/lsw_bioblend.log 2>&1
        sleep 30
    done
) &

exec shiny-server
