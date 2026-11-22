# Public web page (developer notes)

The pipeline can also be run from a public web page, without any account. The design is the same as [SHAMAN](https://github.com/aghozlane/shaman):

<pre>
browser -> web page (web/app/app.R, Shiny)  -- writes -->  todo/KEY.json + jobs/KEY/{sample/*, params.json}
                                                                  |
           web/daemon/lsw_bioblend.py (same container, Galaxy service account, BioBlend)
                                                                  |
           galaxy.pasteur.fr: lab_strain_watch tool (galaxy/lab_strain_watch.xml)  -- Slurm -->  Maestro
</pre>

- The sample files (or the run accession) are sent to the Galaxy tool ([galaxy/README.md](../galaxy/README.md)) as its sample input, and the values of the web form as its JSON parameter file.
- *web/app/app.R*: the web page. The form is built from *nextflow.config* at start-up, one field per parameter (except those of *galaxy/locked_params.txt*), with the comment of the parameter as help. Users get a key to follow the job and download the results: strain identity report, coverage report and zip of the result folder ("Check my job" tab, or link in the e-mail).
- *web/daemon/lsw_bioblend.py*: sends each job to Galaxy, waits, downloads the results, sends an e-mail and removes old jobs. Settings are environment variables (see the top of the file). The Galaxy API key must never be written in the code.
- *web/Dockerfile*: image of the web page and the daemon. LSW_MAX_UPLOAD_MB (1000 by default) limits the total size of the sample files of a submission: sequencing runs are large, and the run accession field of the page avoids their upload (the Galaxy job then needs an internet access, see [galaxy/README.md](../galaxy/README.md)).
- *.gitlab-ci.yml* (root of the repository): builds this image and pushes it into the GitLab registry of a gitlab.pasteur.fr mirror of the repository, from which the Pasteur hosting deploys it.

```
docker build -f web/Dockerfile -t lab_strain_watch_web .
docker run -p 3838:3838 -e GALAXY_API_KEY=... -e LSW_PUBLIC_URL=https://... -v lsw_work:/srv/lsw_work lab_strain_watch_web
```
