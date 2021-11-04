FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG SHASTAVER="0.8.0"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Shasta"
LABEL software.version="9bab4d8"
LABEL description="De novo assembly from Oxford Nanopore reads."
LABEL website="https://github.com/chanzuckerberg/shasta"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y curl ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install shasta binary
RUN cd /usr/bin/ && curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/0.8.0/shasta-Linux-0.8.0 && \
    chmod ugo+x shasta-Linux-0.8.0 && ln -s shasta-Linux-0.8.0 shasta

CMD [ 'shasta' ]
