FROM nvidia/cuda:11.7.0-base-ubuntu22.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG GUPPYVER=6.3.2

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install guppy
RUN apt-get update && \
    apt-get install -y wget && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    tar -zxf ont-guppy_${GUPPYVER}_linux64.tar.gz && rm ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    mv ont-guppy /opt/ && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/ont-guppy/bin/:${PATH}"

CMD [ 'guppy_basecaller' ]