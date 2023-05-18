FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build
ARG SAMTOOLSVER=1.17
ARG DORADOVER=0.2.2
ARG MINIMAP2VER=2.24

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="Samtools"
LABEL software.version="1.15.1"
LABEL description="Tools (written in C using htslib) for manipulating next-generation sequencing data"
LABEL website="https://github.com/samtools/samtools"
LABEL license="https://github.com/samtools/samtools/blob/develop/LICENSE"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 wget \
 make \
 bc \
 perl \
 bzip2 \
 gnuplot \
 libqt5core5a \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
    rm samtools-${SAMTOOLSVER}.tar.bz2 && \
    cd samtools-${SAMTOOLSVER} && \
    ./configure && \
    make && \
    make install

# get plot-bamstats to work
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADOVER}-linux-x64.tar.gz && \
    tar -xf dorado-${DORADOVER}-linux-x64.tar.gz && \
    rm dorado-${DORADOVER}-linux-x64.tar.gz && mv dorado* /opt/dorado && \
    cd /opt/dorado/bin && ./dorado download --model all

ENV PATH=/opt/dorado/bin/:$PATH

# install minimap2 binary
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2VER}/minimap2-${MINIMAP2VER}_x64-linux.tar.bz2 && \
    tar -jxvf minimap2-${MINIMAP2VER}_x64-linux.tar.bz2

ENV PATH="${PATH}:/minimap2-${MINIMAP2VER}_x64-linux"

# set perl locale settings
ENV LC_ALL=C
