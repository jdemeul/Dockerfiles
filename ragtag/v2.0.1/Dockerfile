FROM zeunas/minimap2:2.22

# for easy upgrade later. ARG variables only persist during image build
ARG RAGTAGVER=2.0.1

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="RagTag"
LABEL software.version="2.0.1"
LABEL description="Tools for fast and flexible genome assembly scaffolding and improvement"
LABEL website="https://github.com/malonge/RagTag"
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
    pip install --no-cache-dir RagTag==${RAGTAGVER}

ENV PATH="$HOME/.local/bin:$PATH"

CMD [ "ragtag.py" ]