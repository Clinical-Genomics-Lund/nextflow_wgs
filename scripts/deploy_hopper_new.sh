#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#DIR="/fs1/jakob/deploy/cmd_wgs"


# FIXME: How to run either from hopper or remotely? A prompt maybe
DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/pipelines/wgs_germline38"

current_computer=$(hostname)

# if [[ "${current_computer}" != "rs-fe1" ]]; then
#     DEST="${DEST_HOST}:${PIPELINE_DEST}"
#     echo "You are not currently on host rs-fe1, assuming you are deploying remotely"
#     echo "Destination: ${DEST}"
# else
#     DEST="${PIPELINE_DEST}"
#     echo "You are on host rs-fe1, assuming you are deploying locally"
#     echo "Destination: ${DEST}"
# fi

# rs-fe1

current_branch=$(git branch | grep "^*" | sed "s/^* //")
if [[ "${current_branch}" != "master" ]]; then
    echo "You are not on the master branch. Current branch: ${current_branch}"
    echo "Do you want to deploy anyway? (y/n)"
    read -r deploy_despite_not_master
    if ! [[ "${deploy_despite_not_master}" =~ ^[yY]$ ]]; then
        echo "Aborting"
        exit 0
    fi
fi

echo "Confirm, there are no running jobs for this pipeline? (y/n)"
read -r no_running

if ! [[ "${no_running}" =~ ^[yY]$ ]]; then
    echo "Aborting"
    exit 0
fi

# echo "Retrieving git.hash ..."
# if [[ "${current_computer}" != "rs-fe1" ]]; then
echo "Retrieving current git.hash from ${DEST_HOST}:$PIPELINE_DEST/git.hash ..."
current_hash=$(ssh "${DEST_HOST}" "cat $PIPELINE_DEST/git.hash")
# else
#     current_hash=$(cat $PIPELINE_DEST/git.hash)
# fi

echo "Showing diff to ${current_hash}"
git diff ${current_hash} --stat

echo "Do you want to view the full diff? (y/n)"
read -r view_full

if [[ "${view_full}" =~ ^[Yy]$ ]]; then
    git diff ${current_hash}
fi

echo "Do you want to proceed with deploying? (y/n)"
read -r response

if [[ "${response}" =~ ^[Yy]$ ]]; then
    # Copy pipeline script
    cp -v $DIR/main.nf $DEST

    # Copy configuration file
    cp -v $DIR/configs/nextflow.hopper.config $DEST/nextflow.config

    # Copy other files
    cp -r -v $DIR/bin $DEST

    #git rev-parse HEAD > git.hash
    cp -v $DIR/git.hash $DEST

else
    echo "Deploy aborted"
    exit 0
fi
