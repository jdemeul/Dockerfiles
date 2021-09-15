FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG HTSLIBVER=1.13
ARG SHAPEITVER=4.2.1

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="SHAPEIT"
LABEL software.version="4.2.2"
LABEL description="SHAPEIT is a fast and accurate method for estimation of haplotypes (aka phasing) from genotype or sequencing data."
LABEL website="https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html"
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
 libboost-all-dev \
 g++ \
 wget \
 make \
 bzip2 \
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

# install SHAPEIT
RUN wget https://github.com/odelaneau/shapeit4/archive/refs/tags/v${SHAPEITVER}.tar.gz && \
 tar -xvf v${SHAPEITVER}.tar.gz && \
 rm v${SHAPEITVER}.tar.gz && \
 cd shapeit4-${SHAPEITVER} && \
 sed -i 's/$(HOME)\/Tools\/htslib-1.11/\/htslib-1.13/g' makefile && \
 make && ln -s /shapeit4-${SHAPEITVER}/bin/shapeit4* /usr/bin/shapeit4

CMD [ "shapeit4" ]
