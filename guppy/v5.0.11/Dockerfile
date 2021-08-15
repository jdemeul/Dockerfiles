FROM nvidia/cuda:11.2.1-base-ubuntu20.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG GUPPYVER=5.0.11

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install guppy
RUN apt-get update && \
    apt-get install -y wget && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    tar -zxf ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    mv ont-guppy /opt/ && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/ont-guppy/bin/:${PATH}"

CMD [ 'guppy_basecaller' ]