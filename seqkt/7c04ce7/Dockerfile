FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG SEQTKVER="7c04ce7"

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="vcf2diploid"
LABEL software.version="d4517c5"
LABEL description="CrossStitch creates personalized reference-quality diploid genomes without de novo assembly."
LABEL website="https://github.com/schatzlab/crossstitch"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install -y \
 git \
 make \
 gcc \
 zlib1g-dev \
 ca-certificates \
 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install HapCUT2
RUN git clone https://github.com/lh3/seqtk && \
    cd seqtk && git checkout -b ${SEQTKVER} && \
    make && mv seqtk /usr/local/bin/. && rm -rf /seqtk

CMD [ "seqtk" ]
