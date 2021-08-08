FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG FASTQCVER=0.11.9

LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="FASTQC"
LABEL software.version="0.11.9"
LABEL about.summary="A quality control tool for high throughput sequence data."
LABEL about.home="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
LABEL about.documentation="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/"
LABEL license="https://www.gnu.org/copyleft/gpl.html"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v${FASTQCVER}.zip

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre perl unzip \
    && cd /usr/local/ \
    && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${ZIP} \
    && unzip ${ZIP} \
    && rm -rf ${ZIP} \
    && chmod 755 FastQC/fastqc \
    && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set perl locale settings
ENV LC_ALL=C

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/fastqc"]