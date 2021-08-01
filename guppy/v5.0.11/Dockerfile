FROM nvidia/cuda:11.4.0-base-ubuntu20.04

LABEL base.image="ubuntu:xenial"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG GUPPYVER=v5.0.11

# install fastp
RUN apt-get update && \
    apt-get install -y wget lsb-release && \
    export PLATFORM=$(lsb_release -cs) && \
    wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - && \
    echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list && \
    apt-get update && \
    apt install -y ont-guppy && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*