FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG SHASTAVER="0.11.1"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Shasta"
LABEL software.version="9a5cfdd"
LABEL description="De novo assembly from Oxford Nanopore reads."
LABEL website="https://github.com/paoloshasta/shasta"
LABEL maintainer="Alexandra Pancikova"
LABEL maintainer.email="alexandra.pancikova@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y curl git ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install shasta binary
RUN cd /usr/bin/ && curl -O -L https://github.com/paoloshasta/shasta/releases/download/${SHASTAVER}/shasta-Linux-${SHASTAVER} && \
    chmod ugo+x shasta-Linux-${SHASTAVER} && ln -s shasta-Linux-${SHASTAVER} shasta


CMD [ "shasta" ]