#!/bin/bash
set -e
WORKDIR=lfa_lab
(cd $WORKDIR; git pull)
(cd $WORKDIR; cmake . ; make -j4 sphinx-doc)
# copy files and delete old ones
# keep all configuration files and files we created in this repository
rsync -rv --delete \
  --exclude='update.sh' \
  --exclude='releases' \
  --exclude='lfa_lab' \
  --exclude='.*' \
  $WORKDIR/doc/html/ .

