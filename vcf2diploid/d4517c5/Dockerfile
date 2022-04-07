FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG VCF2DIPLOIDVER="d4517c5"

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="vcf2diploid"
LABEL software.version="d4517c5"
LABEL description="CrossStitch creates personalized reference-quality diploid genomes without de novo assembly."
LABEL website="https://github.com/schatzlab/crossstitch"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 git \
 make \
 ca-certificates \
 default-jdk \
 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install HapCUT2
RUN git clone https://github.com/abyzovlab/vcf2diploid && \
    cd vcf2diploid && \
    make

ENTRYPOINT ["java","-Xmx31G","-jar","/vcf2diploid/vcf2diploid.jar"]
