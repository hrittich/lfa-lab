#!/bin/bash
set -e
WORKDIR=$(mktemp -d)
git clone https://github.com/hrittich/lfa-lab $WORKDIR
(cd $WORKDIR; cmake . ; make -j4 sphinx-doc)
rsync -rv --delete \
  --exclude='update.sh' \
  --exclude='releases' \
  --exclude='.*' \
  $WORKDIR/doc/html/ .

rm -r $WORKDIR

