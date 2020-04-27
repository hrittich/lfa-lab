#!/bin/bash
set -xe

SHELL_ON_ERROR=${SHELL_ON_ERROR:-yes}

git bundle create repository.git HEAD

TESTS="python2-out-of-source python3-in-source python3-pip python3-venv"
for TEST in $TESTS; do
  CID=$(docker create -ti --rm -e SHELL_ON_ERROR=$SHELL_ON_ERROR lfa-lab-testing $TEST)
  docker cp repository.git $CID:/home/tester/repository.git
  docker start -i $CID
done

