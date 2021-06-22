#!/bin/bash
set -x

case "$1" in
  shell)
    exec bash
    ;;
  *)
    git clone /source ~tester/source
    chown -R tester:tester ~tester/source
    sudo -i -u tester bash ~tester/source/docker/run-test-setting.sh "$@"
  ;;
esac

