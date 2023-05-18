FROM nvidia/cuda:11.2.2-cudnn8-runtime-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Jonas Demeulemeester"
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG MEDAKAVER=1.8.0

RUN \
    apt update && \
    apt install -yq --no-install-recommends \
        ca-certificates build-essential cmake curl wget git git-lfs \
        zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev \
        libssl-dev libffi-dev \
        libreadline8 libreadline-dev sqlite3 libsqlite3-dev file \
        python3-all-dev python3-venv python3-pip python3-setuptools && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 10 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN git lfs install

RUN git clone --branch v${MEDAKAVER} --depth 1 https://github.com/nanoporetech/medaka.git && \
    cd medaka && \
    make install_root && \
    cd / && \
    rm -rf /medaka

CMD [ 'medaka' ]