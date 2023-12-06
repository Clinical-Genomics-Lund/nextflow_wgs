# Deploy

Only a subset of the repository is needed to deploy.

Key files are:

* `main.nf`: The main Nextflow script
* `nextflow.config`: Configurations (location of annotations etc.)
* `git.hash`: Hash for latest commit to display in log file
* `bin`-folder: Collection of utility scripts

Steps to deploy:

* Prior to deployment, make sure there are no currently running jobs using the workflow
* Deploy by running the `deploy_hopper_dev.sh` script (check that `DEST_HOST` and `PIPELINE_DEST` are what you would expect)
* Verify that deployment has been successful by checking that the deployed `git.hash` contains the expected ID
