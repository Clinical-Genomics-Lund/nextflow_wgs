DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

LATEST_CONTAINER_BUILD="$( ls -t $DIR/container/wgs_*.sif |head -n1)"
CONTAINER_BASENAME=${LATEST_CONTAINER_BUILD##*/}
PIPELINE_DEST="/media/hopper/pipelines/wgs_germline"
CONTAINER_DEST=/media/hopper/resources/containers/$CONTAINER_BASENAME


# Deploy container if it isn't already deployed
if test -f "$CONTAINER_DEST"; then
    echo "Latest container already deployed, skipping!"
else
    echo "Deploying container"
    cp $LATEST_CONTAINER_BUILD $CONTAINER_DEST
    # TODO: Replace "active" container symlink on hopper!
fi


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp $DIR/shards.csv $PIPELINE_DEST
cp $DIR/test_run.sh $PIPELINE_DEST
cp $DIR/params.json $PIPELINE_DEST
