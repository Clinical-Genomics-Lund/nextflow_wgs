DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

LATEST_CONTAINER_BUILD="$( ls -t container/wgs_*.sif |head -n1)"

PIPELINE_DEST="/media/hopper/pipelines/wgs_germline"

# Copy container
#cp $DIR/$LATEST_CONTAINER_BUILD /media/hopper/resources/containers/


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp $DIR/shards.csv $PIPELINE_DEST
cp $DIR/test_run.sh $PIPELINE_DEST
cp $DIR/params.json $PIPELINE_DEST
