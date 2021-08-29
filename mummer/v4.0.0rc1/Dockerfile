FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG MUMMERVER="4.0.0rc1"

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="MUMmer"
LABEL software.version="e49527b"
LABEL description="MUMmer is a versatile alignment tool for DNA and protein sequences."
LABEL website="https://github.com/mummer4/mummer"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 build-essential \
 wget \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install MUMmer
RUN wget https://github.com/mummer4/mummer/releases/download/v${MUMMERVER}/mummer-${MUMMERVER}.tar.gz && \
    tar -zxf mummer-${MUMMERVER}.tar.gz && \
    rm mummer-${MUMMERVER}.tar.gz && \
    cd mummer-${MUMMERVER} && \
    ./configure --prefix=/usr/ && \
    make && make install

# set perl locale settings
CMD [ "nucmer" ]
