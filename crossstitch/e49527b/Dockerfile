FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG CROSSSTITCHVER="e49527b"
ARG HAPCUTVER=1.3.3
ARG HTSLIBVER=1.13

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="CrossStitch"
LABEL software.version="e49527b"
LABEL description="CrossStitch creates personalized reference-quality diploid genomes without de novo assembly."
LABEL website="https://github.com/schatzlab/crossstitch"
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
 git \
 g++ \
 wget \
 make \
 perl \
 bzip2 \
 gnuplot \
 default-jre \
 default-jdk \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-${HTSLIBVER}.tar.bz2 && \
 tar -xjf htslib-${HTSLIBVER}.tar.bz2 && \
 rm htslib-${HTSLIBVER}.tar.bz2 && \
 cd htslib-${HTSLIBVER} && \
 ./configure --prefix=/usr/ && \
 make && \
 make install

# install HapCUT2
RUN git clone --depth 1 --branch "v.${HAPCUTVER}" https://github.com/vibansal/HapCUT2 && \
    cd HapCUT2 && \
    make

# install CrossStitch
RUN git clone --recursive https://github.com/schatzlab/crossstitch.git && \
    cd crossstitch && \
    git reset --hard ${CROSSSTITCHVER} && \
    cd vcf2diploid && \
    make && \
    cd ../Iris && \
    ./build.sh && \
    mv /HapCUT2/build/HAPCUT2 /usr/bin && mv /HapCUT2/build/extractHAIRS /usr/bin

# set perl locale settings
CMD [ "/crossstitch/src/crossstitch.sh" ]
