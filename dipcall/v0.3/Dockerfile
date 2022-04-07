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
RUN apt-get update && apt-get install -y wget perl make && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install minimap2 binary
RUN wget https://github.com/lh3/dipcall/releases/download/v${DIPCALLVER}/dipcall-${DIPCALLVER}_x64-linux.tar.bz2 && \
    tar -jxf dipcall-${DIPCALLVER}_x64-linux.tar.bz2

RUN cd dipcall.kit && wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed

# set perl locale settings
ENV LC_ALL=C

ENV PATH="${PATH}:/dipcall.kit/"
CMD [ 'run-dipcall' ]