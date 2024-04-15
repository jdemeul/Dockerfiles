FROM continuumio/miniconda3:4.12.0
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip ca-certificates && \
    apt-get install --no-install-recommends -y lsb-release gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN conda install python=3 numpy scipy && \
    conda install -c bioconda minimap2 samtools tabix spoa longshot bcftools qualimap

RUN cd /opt && git clone -b v1.1 https://github.com/gaolabtools/scNanoGPS.git && \
    cd scNanoGPS && pip3 install -r requirements.txt

RUN pip3 install liqa

# RUN wget 

ENV PATH=/opt/scNanoGPS:$PATH
