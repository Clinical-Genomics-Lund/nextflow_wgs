DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/pipelines/wgs_germline_oncotmp"

# Copy pipeline script
scp $DIR/main.nf $DEST_HOST:$PIPELINE_DEST

# Copy configuration file
scp $DIR/configs/nextflow.hopper.config $DEST_HOST:$PIPELINE_DEST/nextflow.config

# Copy other files
scp $DIR/shards_38.csv $DEST_HOST:$PIPELINE_DEST
scp -r $DIR/bin $DEST_HOST:$PIPELINE_DEST
