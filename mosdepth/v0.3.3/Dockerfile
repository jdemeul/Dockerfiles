FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG MOSDEPTHVER=0.3.3

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="mosdepth"
LABEL software.version="0.3.3"
LABEL description="A tool for quick quality assessment of cram and bam files, intended for long read sequencing."
LABEL website="https://github.com/brentp/mosdepth"
LABEL license="https://github.com/brentp/mosdepth/blob/master/LICENSE"
LABEL maintainer="Alexandra Pancikova"
LABEL maintainer.email="alexandra.pancikova@kuleuven.be"

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
 zip \
 unzip \ 
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTHVER}/mosdepth && \
    chmod +x mosdepth && \
    mv mosdepth /usr/local/bin/mosdepth


ENV PATH="${PATH}:/mosdepth/"



# set perl locale settings
ENV LC_ALL=C