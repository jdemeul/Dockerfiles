FROM nvidia/cuda:10.1-base-ubuntu18.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG MEDAKAVER=1.4.3
ARG HTSLIBVER=1.13
ARG SAMTOOLSVER=1.13
ARG MINIMAPVER=2.21

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 libffi-dev \
 libncurses5-dev \
 gcc \
 g++ \
 wget \
 make \
 cmake \
 ca-certificates \
 python3 \
 python3-all-dev \
 python3-pip \
 curl \
 bzip2 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-${HTSLIBVER}.tar.bz2 && \
 tar -xjf htslib-${HTSLIBVER}.tar.bz2 && \
 rm htslib-${HTSLIBVER}.tar.bz2 && \
 cd htslib-${HTSLIBVER} && \
 ./configure --prefix=/usr/ && \
 make && \
 make install

# install minimap2 binary
RUN curl -L https://github.com/lh3/minimap2/releases/download/v${MINIMAPVER}/minimap2-${MINIMAPVER}_x64-linux.tar.bz2 | tar -jxvf - && \
  chown root:root minimap2-${MINIMAPVER}_x64-linux && mv minimap2-${MINIMAPVER}_x64-linux /usr/bin/minimap2

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir medaka==${MEDAKAVER}