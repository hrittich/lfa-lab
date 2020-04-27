#!/bin/bash
set -e

run_test_suite () {
  echo "### running test suite: $1 ###"
  set -xe
  case "$1" in
  shell)
    bash

    ;;
  python2-out-of-source)
    mkdir build
    (cd build &&
     cmake -DPYTHON_EXECUTABLE=`which python2` -DWITH_TESTS=ON \
       -DUSER_INSTALL=ON ../ )
    make -C build
    make -C build check
    make -C build install
    (cd demo; python2 -mlfa_lab)

    ;;
  python3-in-source)

    cmake -DPYTHON_EXECUTABLE=`which python3` -DWITH_TESTS=ON \
      -DUSER_INSTALL=ON .
    make
    make check
    make install
    (cd demo; python3 -mlfa_lab)

    ;;
  python3-pip)

    pip3 install --no-deps -v --user .
    (cd demo && python3 -mlfa_lab)

    ;;
  python3-venv)

    python3 -mvenv myenv
    (source myenv/bin/activate; pip install -v .)
    (source myenv/bin/activate; cd demo; python -mlfa_lab)

    ;;
  *)
    echo "Invalid test suite name"
    exit 1
  esac
}

git clone repository.git build
cd build
git checkout -b testing

# Weird POSIX. See: https://unix.stackexchange.com/questions/65532/
set +e
(run_test_suite "$1")
EXIT_STATUS=$?
set -e

if [ $EXIT_STATUS -ne 0 ]; then
  [ "$SHELL_ON_ERROR" == yes ] && bash
fi

