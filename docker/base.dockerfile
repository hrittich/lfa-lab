FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

# Install dev tool needed for LFA Lab
RUN \
  apt-get update && apt-get install -q -y \
    cmake \
    g++ \
    git \
    libeigen3-dev \
    libgtest-dev \
    liblapack-dev \
    python3-dev \
    python3-matplotlib \
    python3-numpy-dev \
    python3-pip \
    python3-setuptools \
    python3-six \
    python3-venv \
    python-dev \
    python-matplotlib \
    python-numpy-dev \
    python-pip \
    python-setuptools \
    python-six \
    sudo \
    swig \
    vim \
    wget

# get version 3.12.4 of cmake
WORKDIR /root
RUN \
  wget --progress=dot:mega https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz 2>&1 && \
  echo "486edd6710b5250946b4b199406ccbf8f567ef0e23cfe38f7938b8c78a2ffa5f cmake-3.12.4-Linux-x86_64.tar.gz" \
    | sha256sum -c - && \
  tar -C /usr/local --strip-components=1 -zxf cmake-3.12.4-Linux-x86_64.tar.gz && \
  cmake -version && \
  rm cmake-3.12.4-Linux-x86_64.tar.gz

CMD ["/bin/bash"]
