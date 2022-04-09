FROM continuumio/miniconda3:4.10.3
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN cd /opt && \
    wget -O minimap2.tar.gz https://github.com/lh3/minimap2/archive/v2.24.tar.gz && \
    mkdir minimap2 && \
    tar -zxvf minimap2.tar.gz -C minimap2 --strip-components 1 && \
    rm minimap2.tar.gz && \
    cd minimap2 && \
    make

ENV PATH=/opt/minimap2/:$PATH

RUN cd /opt && \
    wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    mkdir samtools && \
    tar -xvjf samtools.tar.bz2 -C samtools --strip-components 1 && \
    rm samtools.tar.bz2 && \
    cd samtools && \
    ./configure && \
    make && \
    make install

RUN conda install -y -c defaults -c conda-forge -c bioconda pysam

COPY haplotaggedbam_to_fastq.py /usr/local/bin/
COPY retag_bam.py /usr/local/bin/

RUN chmod +x /usr/local/bin/haplotaggedbam_to_fastq.py /usr/local/bin/retag_bam.py
