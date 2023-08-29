FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build
ARG HIFIASMVER=0.19.5

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="Samtools"
LABEL software.version="1.18"
LABEL description="Tools (written in C using htslib) for manipulating next-generation sequencing data"
LABEL website="https://github.com/samtools/samtools"
LABEL license="https://github.com/samtools/samtools/blob/develop/LICENSE"
LABEL maintainer="Alexandra Pancikova"
LABEL maintainer.email="alexandra.pancikova@kuleuven.be"


ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install -y autoconf --no-install-recommends -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 g++ \
 wget \
 make \
 bc \
 bzip2 \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${HIFIASMVER}.tar.gz && \
    tar -xzf ${HIFIASMVER}.tar.gz && \
    rm ${HIFIASMVER}.tar.gz && \
    cd hifiasm-${HIFIASMVER} && \
    make

RUN mv hifiasm-${HIFIASMVER}/hifiasm /usr/local/bin/ && \
    rm -rf hifiasm-${HIFIASMVER}

