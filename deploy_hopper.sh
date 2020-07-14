DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

LATEST_CONTAINER_BUILD="$( ls -t $DIR/container/wgs_*.sif |head -n1)"
CONTAINER_BASENAME=${LATEST_CONTAINER_BUILD##*/}
PIPELINE_DEST="/fs1/pipelines/wgs_germline38"
CONTAINER_DEST=/fs1/resources/containers/$CONTAINER_BASENAME
DEST_HOST="rs-fs1.lunarc.lu.se"

read -p "Are you sure you want to deploy to production? " -n 1 -r

echo
if [[ $REPLY =~ ^[Yy]$ ]]
then

    # Deploy container if it isn't already deployed
    if test -f "$CONTAINER_DEST"; then
	echo "Latest container already deployed, skipping!"
    else
	echo "Deploying container"
	scp $LATEST_CONTAINER_BUILD $DEST_HOST:$CONTAINER_DEST
	# TODO: Replace "active" container symlink on hopper!
    fi


    # Copy pipeline script
    scp $DIR/main.nf $DEST_HOST:$PIPELINE_DEST

    # Copy configuration file
    scp $DIR/configs/nextflow.hopper.config $DEST_HOST:$PIPELINE_DEST/nextflow.config

    # Copy other files
    scp $DIR/shards.csv $DEST_HOST:$PIPELINE_DEST
    scp -r $DIR/bin $DEST_HOST:$PIPELINE_DEST
fi 
