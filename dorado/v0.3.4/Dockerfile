FROM nvidia/cuda:12.2.0-base-ubuntu22.04

# for easy upgrade later. ARG variables only persist during image build
ARG DORADOVER=0.3.4

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
    cd /opt/dorado/bin && ./dorado download --model all

ENV PATH=/opt/dorado/bin/:$PATH

# set perl locale settings
ENV LC_ALL=C
