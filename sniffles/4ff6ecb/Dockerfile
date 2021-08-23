FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG SNIFFLESVER="4ff6ecb"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Sniffles"
LABEL software.version="4ff6ecb"
LABEL description="Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore)."
LABEL website="https://github.com/fritzsedlazeck/Sniffles"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 libffi-dev \
 libncurses5-dev \
 git \
 file \
 gcc \
 g++ \
 wget \
 make \
 cmake \
 curl \
 bzip2 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install sniffles binary
RUN git clone https://github.com/fritzsedlazeck/Sniffles && cd Sniffles && mkdir -p build && cd build && cmake .. && make && \
    mv /Sniffles/bin/sniffles*/sniffles* /usr/bin/. && rm -rf /Sniffles

CMD [ 'sniffles' ]