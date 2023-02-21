FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build
ARG CRAMINOVER=0.9.6

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="cramino"
LABEL software.version="0.9.6"
LABEL description="A tool for quick quality assessment of cram and bam files, intended for long read sequencing."
LABEL website="https://github.com/wdecoster/cramino"
LABEL license="https://github.com/wdecoster/cramino/blob/master/LICENSE"
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

RUN wget https://github.com/wdecoster/cramino/releases/download/v${CRAMINOVER}/cramino-linux.zip && \
    unzip cramino-linux.zip && \
    chmod +x cramino && \
    mv cramino /usr/local/bin/cramino


ENV PATH="${PATH}:/cramino/"



# set perl locale settings
ENV LC_ALL=C



