FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG HTSLIBVER=1.14
ARG LRAVER=1.3.3

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="LRA"
LABEL software.version="1.3.3"
LABEL description="lra is a sequence alignment program that aligns long reads."
LABEL website="https://github.com/ChaissonLab/LRA"
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
 g++ \
 wget \
 make \
 perl \
 bzip2 \
 gnuplot \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIBVER}/htslib-${HTSLIBVER}.tar.bz2 && \
 tar -xjf htslib-${HTSLIBVER}.tar.bz2 && \
 rm htslib-${HTSLIBVER}.tar.bz2 && \
 cd htslib-${HTSLIBVER} && \
 ./configure --prefix=/usr/ && \
 make && \
 make install

# install LRA
RUN wget https://github.com/ChaissonLab/LRA/archive/refs/tags/v${LRAVER}.tar.gz && \
 tar -xvf v${LRAVER}.tar.gz && \
 rm v${LRAVER}.tar.gz && \
 cd LRA-${LRAVER} && \
 make && \
 mv lra /usr/bin && cd .. && rm -rf LRA-${LRAVER}

# set perl locale settings
ENV LC_ALL=C