FROM ubuntu:24.04

ENV PATH="/fgbio/miniconda3/bin:${PATH}"
ARG PATH="/fgbio/miniconda3/bin:${PATH}"

RUN apt-get update

RUN apt-get install -y wget bzip2 && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir -p /fgbio/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /fgbio/miniconda3 -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | \
#     tar -xvj bin/micromamba && \
#     bin/micromamba shell init -s bash -r micromamba

# RUN micromamba config append channels conda-forge && \
#     micromamba config append channels nodefaults && \
#     micromamba config set auto_activate_base true

RUN conda config --add channels bioconda
RUN conda install -y fgbio