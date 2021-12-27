FROM ubuntu:focal

ARG SAMTOOLSVER=1.14
ARG DIPDIFFVER=c773070

# for easy upgrade later. ARG variables only persist during image build
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# update and install dependencies
# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 g++ \
 git \
 python3-pip \
 wget \
 make \
 perl \
 bzip2 \
 tzdata \
 gnuplot \
 libqt5core5a \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install && \
 cd ..

RUN git clone https://github.com/fenderglass/dipdiff && \
    cd dipdiff && \
    git checkout ${DIPDIFFVER} && \
    git submodule update --init && \
    make && \
    python3 -m pip install --upgrade pip && \
    pip install -r requirements.txt

ENV PATH "/dipdiff/:/dipdiff/submodules/minimap2/:/dipdiff/submodules/svim-asm/src/svim_asm/:${PATH}"
ENV PYTHONUNBUFFERED "1"
ENV MPLCONFIGDIR "/tmp"