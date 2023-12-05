* Key files to deploy
    * `main.nf`
    * `configs/nextflow.hopper.config` -> `nextflow.config` next to `main.nf`
    * `shards.csv` / `shards_38.csv`
    * `git.hash`
    * `bin`-folder
* Update the config
* Place these files on hopper
* Files can be deployed by running `deploy_hopper.sh` / `deploy_hopper_dev.sh`

* Running using the Hopper Perl script
    * Specify "nolog"
    * Specify output folder
