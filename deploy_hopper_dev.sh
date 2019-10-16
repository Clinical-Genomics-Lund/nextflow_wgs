DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="/media/hopper/pipelines/wgs_germline_dev/"


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp $DIR/shards.csv $PIPELINE_DEST
cp -r $DIR/bin $PIPELINE_DEST
cp -r $DIR/rank_models /fs1/resources/scout/
