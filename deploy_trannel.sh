DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="/trannel/pipelines/onco-solid"


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.trannel.config $PIPELINE_DEST/nextflow.config
cp $DIR/shards_38.csv $PIPELINE_DEST
# Update bin-folder from github
#. bin-update.sh
# Copy other files
cp -r $DIR/bin $PIPELINE_DEST

