FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG WHATSHAPVER=1.1

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="WhatsHap"
LABEL software.version="1.1"
LABEL description="lra is a sequence alignment program that aligns long reads."
LABEL website="https://whatshap.readthedocs.io/en/latest/index.html"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 build-essential \
 python3-dev \
 python3-pip && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir whatshap==${WHATSHAPVER}

ENV PATH="$HOME/.local/bin:$PATH"

CMD [ "whatshap" ]