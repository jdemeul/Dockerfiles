FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG PIGZVER=2.6
ARG ZLIBNGVER=2.0.5

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Samtools"
LABEL software.version="1.13"
LABEL description="Tools (written in C using htslib) for manipulating next-generation sequencing data"
LABEL website="https://github.com/samtools/samtools"
LABEL license="https://github.com/samtools/samtools/blob/develop/LICENSE"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 gcc \
 g++ \
 wget \
 make \
 cmake \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install pigz zlib-ng
RUN wget https://github.com/zlib-ng/zlib-ng/archive/refs/tags/${ZLIBNGVER}.tar.gz && \
 tar -zxf ${ZLIBNGVER}.tar.gz && \
 cd zlib-ng-${ZLIBNGVER} && \
 cmake . -DZLIB_COMPAT=ON && cmake --build . --target install && \
 ldconfig
#  ./configure && \
# #  make && make install

RUN wget https://github.com/madler/pigz/archive/refs/tags/v${PIGZVER}.tar.gz && \
 tar -xzf v${PIGZVER}.tar.gz && \
 cd pigz-${PIGZVER} && make && \
 mv pigz /usr/bin/. && mv unpigz /usr/bin/. && \
 cd .. && rm -rf zlib-ng-${ZLIBNGVER} pigz-${PIGZVER}

CMD [ 'pigz' ]