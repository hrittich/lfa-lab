#!/bin/bash
set -e

docker build -f base.dockerfile --tag="lfa-lab-base" .
docker build -f testing.dockerfile --tag="lfa-lab-testing" .

