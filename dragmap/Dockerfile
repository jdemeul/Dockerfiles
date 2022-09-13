FROM continuumio/miniconda3:4.12.0
SHELL ["/bin/bash", "-c"]

# RUN apt-get update && \
#     apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip && \
#     apt-get -qq -y autoremove && \
#     apt-get clean

RUN conda install -y -c defaults -c conda-forge -c bioconda dragmap

