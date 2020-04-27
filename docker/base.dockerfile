FROM ubuntu:bionic

SHELL ["/bin/bash", "-c"]

ENV DEBIAN_FRONTEND=noninteractive

# Install basic development tools
RUN apt-get update
RUN apt-get -y install build-essential

# Install dev tool needed for LFA Lab
RUN \
  apt-get install -y g++ cmake libgtest-dev libeigen3-dev swig liblapack-dev \
      python3-dev python3-numpy-dev python3-matplotlib python3-six \
      python3-pip python3-venv python3-setuptools \
      python-dev python-numpy-dev python-matplotlib python-six \
      python-pip python-setuptools

RUN \
  apt-get install -y git vim sudo wget

# get version 3.12.4 of cmake
WORKDIR /root
RUN \
  wget --progress=dot:mega https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz 2>&1 && \
  echo "486edd6710b5250946b4b199406ccbf8f567ef0e23cfe38f7938b8c78a2ffa5f cmake-3.12.4-Linux-x86_64.tar.gz" \
    | sha256sum -c -
RUN \
  tar -C /usr/local --strip-components=1 -zxf cmake-3.12.4-Linux-x86_64.tar.gz && \
  cmake -version

CMD ["/bin/bash"]

