#!/bin/bash
set -xe
SCRIPTDIR="$(dirname "$0")"

DOCKER=${DOCKER:-docker}
SHELL_ON_ERROR=${SHELL_ON_ERROR:-yes}

ALL_TESTS="python2-out-of-source python3-in-source python3-pip python3-venv"
TESTS="${1:-$ALL_TESTS}"
for TEST in $TESTS; do
  $DOCKER run -ti --rm $DOCKER_FLAGS \
    --env SHELL_ON_ERROR \
    --volume "$SCRIPTDIR/../":/source \
    lfa-lab-testing $TEST
done

