FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG MINIASMVER="ce615d1"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Miniasm"
LABEL software.version="ce615d1"
LABEL description="Miniasm is a very fast OLC-based de novo assembler for noisy long reads."
LABEL website="https://github.com/lh3/miniasm"
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

# install minimap2 binary; make /data
RUN git clone https://github.com/lh3/minimap2 && cd minimap2 && make && mv minimap2 /usr/bin/. && cd .. && \
    git clone https://github.com/lh3/miniasm  && cd miniasm && make && mv miniasm /usr/bin/. && mv minidot /usr/bin/. && cd .. && \
    rm -rf minimap2 miniasm

CMD [ 'miniasm' ]