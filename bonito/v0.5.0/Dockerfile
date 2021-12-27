FROM nvidia/cuda:10.2-cudnn7-devel-ubuntu18.04
# FROM pytorch/pytorch:1.10.0-cuda11.3-cudnn8-runtime

LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG BONITOVER="e4a8e11"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y python3-pip zlib1g-dev git python3-venv && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# # install Bonito
RUN git clone https://github.com/nanoporetech/bonito && \
    cd bonito && \
    git checkout ${BONITOVER}

ENV VIRTUAL_ENV="/bonito/venv3"
RUN cd /bonito && python3 -m venv venv3
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN cd /bonito && \
    pip3 install --no-cache-dir --upgrade pip && \
    pip3 install -r requirements.txt && \
    unset CUDA_VERSION && \
    python3 setup.py develop

# install models
RUN bonito download --models

ENV LC_ALL=C
# LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH" \