FROM continuumio/miniconda3:4.12.0
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip ca-certificates && \
    apt-get install --no-install-recommends -y lsb-release gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN git clone -b v0.2.1 https://github.com/nanoporetech/sockeye.git && \
    cd sockeye && conda env create -f environment.yml

ENV PATH=/opt/conda/envs/sockeye/bin:$PATH
