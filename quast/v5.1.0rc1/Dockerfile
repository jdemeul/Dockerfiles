FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG QUASTVER="5.1.0rc1"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="quast"
LABEL software.version="5.1.0rc1"
LABEL description="QUAST stands for QUality ASsessment Tool. It evaluates genome/metagenome assemblies by computing various metrics"
LABEL website="https://github.com/ablab/quast"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y build-essential pkg-config libfreetype6-dev libpng-dev python3-matplotlib python3-setuptools wget python-is-python3 && \
  apt-get clean && \
  apt-get autoclean && \
  rm -rf /var/lib/apt/lists/*

# for singularity compatibility
ENV LC_ALL=C

# install quast and make /data directory

RUN wget https://github.com/ablab/quast/releases/download/quast_${QUASTVER}/quast-${QUASTVER}.tar.gz && \
    tar -xzf quast-${QUASTVER}.tar.gz && \
    rm -rf quast-${QUASTVER}.tar.gz && \
    cd /quast-${QUASTVER} && \
    /quast-${QUASTVER}/setup.py install_full
