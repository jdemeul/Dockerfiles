FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG SVABAVER=0f60e36

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
 bwa \
 samtools \
 git \
 gcc \
 g++ \
 wget \
 make \
 bzip2 \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/walaj/svaba && \
    cd svaba && git checkout ${SVABAVER} && ./configure && \
    make && make install && \
    mv /svaba/bin/* /usr/local/bin/. && rm -rf /svaba

# install samtools, make /data
# RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
#  tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
#  rm samtools-${SAMTOOLSVER}.tar.bz2 && \
#  cd samtools-${SAMTOOLSVER} && \
#  ./configure && \
#  make && \
#  make install && \
#  mkdir /data

CMD [ 'svaba' ]
