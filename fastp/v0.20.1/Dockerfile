FROM ubuntu:focal

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG FASTPVER=v0.20.1

# install fastp
RUN apt-get update && \
    apt-get install -y git make g++  zlib1g zlib1g-dev
RUN cd /tmp/ && \
    git clone https://github.com/OpenGene/fastp.git && \
    cd fastp && \
    git checkout ${FASTPVER} && \
    make && \
    make install

# clean
RUN rm -rf /tmp/fastp* && \
    apt-get clean && \
    apt-get remove --yes --purge git make g++

CMD fastp