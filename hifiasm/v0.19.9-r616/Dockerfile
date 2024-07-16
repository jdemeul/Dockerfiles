FROM ubuntu:noble

# for easy upgrade later. ARG variables only persist during image build
ARG HIFIASMVER=0.19.9

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="hifiasm"
LABEL software.version="0.19.9"
LABEL description="Hifiasm: a haplotype-resolved assembler for accurate Hifi reads"
LABEL website="https://github.com/chhylp123/hifiasm"
LABEL license="https://github.com/chhylp123/hifiasm/LICENSE"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"


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

