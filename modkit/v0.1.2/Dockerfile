FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build
ARG MODKITVER=0.1.2

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="modkit"
LABEL software.version="0.1.2"
LABEL description="NA"
LABEL website="https://github.com/nanoporetech/modkit/"
LABEL license="https://github.com/nanoporetech/modkit/blob/master/LICENCE.txt"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 wget \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN wget https://github.com/nanoporetech/modkit/releases/download/v${MODKITVER}/modkit_x86_64.tar.gz && \
    tar -xzf modkit_x86_64.tar.gz && \
    rm modkit_x86_64.tar.gz && \
    mkdir -p /opt/modkit/bin && mv dist/modkit /opt/modkit/bin/.

ENV PATH=/opt/modkit/bin/:$PATH

# set perl locale settings
ENV LC_ALL=C
