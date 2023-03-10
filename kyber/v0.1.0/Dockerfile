FROM ubuntu:jammy


ARG KYBERVER=0.1.0

LABEL software="kyber"
LABEL software.version="0.1.0"
LABEL description="A tool for quick quality assessment of cram and bam files, intended for long read sequencing."
LABEL website="https://github.com/wdecoster/kyber"
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

 RUN wget https://github.com/wdecoster/kyber/releases/download/v${KYBERVER}/kyber-linux.zip && \
    unzip kyber-linux.zip && \
    chmod u+x kyber && \
    mv kyber /usr/local/bin/kyber


ENV PATH="${PATH}:/kyber/"



# set perl locale settings
ENV LC_ALL=C