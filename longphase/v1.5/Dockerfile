FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build time
ARG LONGPHASEVER="1.5"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Longphase"
LABEL software.version="f0ba448"
LABEL website="https://github.com/twolinin/longphase"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y curl ca-certificates xz-utils gawk && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install longphase binary
RUN cd /usr/bin/ && curl -O -L https://github.com/twolinin/longphase/releases/download/v${LONGPHASEVER}/longphase_linux-x64.tar.xz && \
    tar -xf longphase_linux-x64.tar.xz && chmod ugo+x longphase_linux-x64 && mv longphase_linux-x64 longphase

CMD [ 'longphase' ]