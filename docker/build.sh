#!/bin/bash
set -e

DOCKER="${DOCKER:-docker}"

$DOCKER build -f base.dockerfile --tag="lfa-lab-base" .
$DOCKER build -f testing.dockerfile --tag="lfa-lab-testing" .

