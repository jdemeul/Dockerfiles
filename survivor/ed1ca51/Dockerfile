FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG SNIFFLESVER="ed1ca51"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="survivor"
LABEL software.version="4ff6ecb"
LABEL description="Toolset for SV simulation, comparison and filtering"
LABEL website="https://github.com/fritzsedlazeck/SURVIVOR"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
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
 make \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install sniffles binary
RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git && cd SURVIVOR/Debug && make && mv SURVIVOR /usr/bin/.

CMD [ 'SURVIVOR' ]