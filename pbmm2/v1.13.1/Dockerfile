FROM ubuntu:24.04

ENV PATH="/pbmm2/miniconda3/bin:${PATH}"
ARG PATH="/pbmm2/miniconda3/bin:${PATH}"

RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir -p /pbmm2/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /pbmm2/miniconda3 -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda config --add channels bioconda
# RUN conda config --add channels conda-forge
RUN conda install -y pbmm2