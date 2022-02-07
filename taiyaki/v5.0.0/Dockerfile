FROM nvidia/cuda:10.2-base-ubuntu18.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG TAIYAKIVER=v5.0.0

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y git python3-pip python3-venv && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install taiyaki
RUN git clone https://github.com/nanoporetech/taiyaki && \
    cd taiyaki && git checkout ${TAIYAKIVER} && cd ..

ENV VIRTUAL_ENV="/taiyaki/venv3"
RUN cd /taiyaki && python3 -m venv venv3
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -e taiyaki/

