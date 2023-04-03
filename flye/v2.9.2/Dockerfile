FROM ubuntu:xenial

# ARG variables only persist during build time. Adding for easy upgrade later
ARG FLYEVER="2.9.2"

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
 python \
 python-pip \
 python-setuptools \
 zlib1g-dev \
 git \
 wget && \
 rm -rf /var/lib/apt/lists/* && apt-get clean && apt-get autoclean

# install flye and update Minimap2
RUN git clone --recursive --branch "${FLYEVER}" https://github.com/fenderglass/Flye.git Flye && \
  cd Flye && \
  python setup.py build && \
  python setup.py install

# just in case for singularity compatibility
ENV LC_ALL=C

CMD [ 'flye' ]