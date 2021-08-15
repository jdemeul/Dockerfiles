FROM nvidia/cuda:11.2.1-base-ubuntu20.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG GUPPYVER=5.0.11
ARG MEGALODONVER=2.3.4

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install guppy
RUN apt-get update && \
    apt-get install -y wget python3-pip && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    tar -zxf ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    mv ont-guppy /opt/ && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install megalodon
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir megalodon==${MEGALODONVER}