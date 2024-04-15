FROM nvidia/cuda:12.3.2-base-ubuntu22.04

# for easy upgrade later. ARG variables only persist during image build
ARG DORADOVER=0.5.3

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 wget \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADOVER}-linux-x64.tar.gz && \
    tar -xf dorado-${DORADOVER}-linux-x64.tar.gz && \
    rm dorado-${DORADOVER}-linux-x64.tar.gz && mv dorado* /opt/dorado && \
    cd /opt/dorado/bin && ./dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0

ENV PATH=/opt/dorado/bin/:$PATH

# set perl locale settings
ENV LC_ALL=C
