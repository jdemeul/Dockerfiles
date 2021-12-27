FROM nvidia/cuda:11.4.2-base-ubuntu20.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG GUPPYVER=6.0.1
ARG MEGALODONVER="48cae7b"
ARG RERIOTAG="abb2d37"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install guppy
RUN apt-get update && \
    apt-get install -y wget git python3-pip && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    tar -zxf ont-guppy_${GUPPYVER}_linux64.tar.gz && \
    mv ont-guppy /opt/ && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install rerio
RUN git clone https://github.com/nanoporetech/rerio && \
    cd rerio && git checkout ${RERIOTAG} && cd .. && \
    rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_v001 && \
    rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 && \
    rerio/download_model.py rerio/basecall_models/res_dna_r941_prom_modbases_5mC_v001 && \
    rerio/download_model.py rerio/basecall_models/res_dna_r103_q20ea_crf_v034 && \
    rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_5hmC_v001

# install megalodon
RUN git clone https://github.com/nanoporetech/megalodon && \
    cd megalodon && git checkout ${MEGALODONVER} && cd .. && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir ont-pyguppy-client-lib==${GUPPYVER} pyparsing==2.4.7 && \
    pip install --no-cache-dir -e megalodon/
