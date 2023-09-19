FROM ubuntu:jammy

# ARG variables only persist during build time. Adding for easy upgrade later
ARG FLYEVER="40b9ea0"

LABEL base.image="ubuntu:xenial"
LABEL dockerfile.version="1"
LABEL software="Flye"
LABEL software.version="2.9.2"
LABEL description="Fast and accurate de novo assembler for single molecule sequencing reads"
LABEL website="https://github.com/fenderglass/Flye"
LABEL license="https://github.com/fenderglass/Flye/blob/flye/LICENSE"
LABEL maintainer1="Jonas Demeulemeester"
LABEL maintainer1.email="Jonas Demeulemeester"

RUN apt-get update && apt-get install -y make \
 gcc \
 g++ \
 python3 \
 python3-pip \
 python3-setuptools \
 zlib1g-dev \
 git \
 wget && \
 rm -rf /var/lib/apt/lists/* && apt-get clean && apt-get autoclean

# install flye and update Minimap2
RUN git clone --recursive https://github.com/fenderglass/Flye.git Flye && \
  cd Flye && \
  git reset --hard ${FLYEVER} && \
  python3 setup.py build && \
  python3 setup.py install

# just in case for singularity compatibility
ENV LC_ALL=C

CMD [ 'flye' ]