FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG MODBAM2BEDVER="v0.5.1"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Modbam2bed"
LABEL software.version="0.5.1"
LABEL website="https://github.com/epi2me-labs/modbam2bed"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

# install deps and cleanup apt garbage
# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 git \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 wget \
 make \
 automake \
 autoconf \
 perl \
 bzip2 \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN git clone -b ${MODBAM2BEDVER} --recursive https://github.com/epi2me-labs/modbam2bed.git && \
    cd modbam2bed && make 

ENV PATH="${PATH}:/modbam2bed"