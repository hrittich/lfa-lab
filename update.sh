#!/bin/bash
set -e
WORKDIR=lfa-lab
(cd $WORKDIR; git pull; cmake . ; make -j4 sphinx-doc)
rsync -rv --delete \
  --exclude='lfa-lab' \
  --exclude='update.sh' \
  --exclude='.*' \
  $WORKDIR/doc/html/ .
