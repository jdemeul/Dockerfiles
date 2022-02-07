FROM nvidia/cuda:11.4.2-base-ubuntu20.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG REMORAVER="569ded0"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y wget git python3-pip python3-venv && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install remora
RUN git clone https://github.com/nanoporetech/remora && \
    cd remora && git checkout ${REMORAVER} && cd ..

ENV VIRTUAL_ENV="/remora/venv3"
RUN cd /remora && python3 -m venv venv3
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install onnxruntime-gpu && \
    pip3 install --no-cache-dir -e remora/
