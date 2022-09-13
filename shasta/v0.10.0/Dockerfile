FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG SHASTAVER="0.10.0"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="Shasta"
LABEL software.version="1be3d4c"
LABEL description="De novo assembly from Oxford Nanopore reads."
LABEL website="https://github.com/chanzuckerberg/shasta"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y curl git ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install shasta binary
RUN cd /usr/bin/ && curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/${SHASTAVER}/shasta-Linux-${SHASTAVER} && \
    chmod ugo+x shasta-Linux-${SHASTAVER} && ln -s shasta-Linux-${SHASTAVER} shasta

# install shasta configs
RUN mkdir shastaconf && cd shastaconf && \
    git init && git remote add origin -f https://github.com/chanzuckerberg/shasta.git && \
    git sparse-checkout set "conf/*" && \
    git pull origin ${SHASTAVER}

CMD [ 'shasta' ]