FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG GFATOOLSVER="0.5"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="gfatools"
LABEL software.version="0.5"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install -y git gcc make zlib1g-dev && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install minimap2 binary; make /data
RUN git clone --branch "v${GFATOOLSVER}" https://github.com/lh3/gfatools && cd gfatools && \
    make && mv gfatools /usr/bin/. && cd .. && rm -rf gfatools

CMD [ 'gfatools' ]