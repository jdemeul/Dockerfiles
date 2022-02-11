FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG DIPCALLVER="0.3"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="dipcall"
LABEL software.version="0.3"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install -y wget perl && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install minimap2 binary; make /data
RUN wget https://github.com/lh3/dipcall/releases/download/v${DIPCALLVER}/dipcall-${DIPCALLVER}_x64-linux.tar.bz2 && \
    tar -jxf dipcall-${DIPCALLVER}_x64-linux.tar.bz2 && ln -s /dipcall.kit/run-dipcall /usr/bin/run-dipcall

CMD [ 'run-dipcall' ]