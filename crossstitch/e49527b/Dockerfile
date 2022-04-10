FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG CROSSSTITCHVER="e49527b"
ARG HAPCUTVER=1.3.3
ARG HTSLIBVER=1.14
ARG BCFTOOLSVER=1.14
ARG SAMTOOLSVER=1.14

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
 libgsl0-dev \
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
 pigz \
 perl \
 libperl-dev \
 default-jre \
 default-jdk \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIBVER}/htslib-${HTSLIBVER}.tar.bz2 && \
    tar -xjf htslib-${HTSLIBVER}.tar.bz2 && \
    rm htslib-${HTSLIBVER}.tar.bz2 && \
    cd htslib-${HTSLIBVER} && \
    ./configure --prefix=/usr/ && \
    make && \
    make install

# install bcftools, make /data
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLSVER}/bcftools-${BCFTOOLSVER}.tar.bz2 && \
    tar -xjf bcftools-${BCFTOOLSVER}.tar.bz2 && \
    rm bcftools-${BCFTOOLSVER}.tar.bz2 && \
    cd bcftools-${BCFTOOLSVER} && \
    ./configure --prefix=/usr --enable-perl-filters --enable-libgsl && \
    make && \
    make install && \
    cd htslib-${BCFTOOLSVER} && make && mv tabix /usr/bin && mv bgzip /usr/bin && \
    cd / && rm -rf bcftools-${BCFTOOLSVER}

# install samtools, make /data
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
    rm samtools-${SAMTOOLSVER}.tar.bz2 && \
    cd samtools-${SAMTOOLSVER} && \
    ./configure && \
    make && \
    make install && \
    mkdir /data

# install HapCUT2
RUN git clone --depth 1 --branch "v.${HAPCUTVER}" https://github.com/vibansal/HapCUT2 && \
    cd HapCUT2 && \
    make && \
    mv /HapCUT2/build/HAPCUT2 /usr/bin && mv /HapCUT2/build/extractHAIRS /usr/bin

# install CrossStitch
RUN git clone --recursive https://github.com/schatzlab/crossstitch.git && \
    cd crossstitch && \
    git reset --hard ${CROSSSTITCHVER} && \
    rm -rf Iris && \
    cd vcf2diploid && \
    make

RUN sed -i 's/100000/250000/g' crossstitch/src/PhaseSVs.java && \
    javac /crossstitch/src/*.java && \
    sed -i 's/javac/#javac/g' /crossstitch/src/crossstitch.sh && \
    sed -i 's/\/work-zfs\/mschatz1\/mkirsche\/github\/HapCUT2\/build\///g' /crossstitch/src/crossstitch.sh

COPY crossstitch_short.sh /crossstitch/src/crossstitch_short.sh
RUN chmod +x /crossstitch/src/crossstitch_short.sh

# set perl locale settings
CMD [ "/crossstitch/src/crossstitch.sh" ]
