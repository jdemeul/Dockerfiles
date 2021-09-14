FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG POMOXISVER=0.3.9

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="pomoxis"
LABEL software.version="0.3.9"
LABEL description="Pomoxis comprises a set of basic bioinformatic tools tailored to nanopore sequencing"
LABEL website="https://github.com/nanoporetech/pomoxis"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
    gcc \
    g++ \
    zlib1g-dev \
    libncurses5-dev \
    python3-all-dev \
    libhdf5-dev \
    libatlas-base-dev \
    libopenblas-base \
    libopenblas-dev \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    make \
    python3-venv \
    python-is-python3 \
    cmake \
    wget \
    bzip2 \
    git \
    ca-certificates && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN git clone --branch v${POMOXISVER} --recursive https://github.com/nanoporetech/pomoxis && \
    cd pomoxis && \
    make install

# add virtualenv to path
ENV PATH="/pomoxis/venv/bin:$PATH"
