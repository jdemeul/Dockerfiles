FROM tensorflow/tensorflow:2.2.2-gpu-py3

LABEL base.image="tensorflow:2.2.2-gpu-py3"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG MEDAKAVER=1.4.3

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 libffi-dev \
 libncurses5-dev \
 git \
 git-lfs \
 file \
 gcc \
 g++ \
 wget \
 make \
 cmake \
 python3-all-dev \
 python3-venv \
 curl \
 bzip2 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install medaka + dependencies
RUN git lfs install

RUN git clone --branch v${MEDAKAVER} --depth 1 https://github.com/nanoporetech/medaka.git && \
    cd medaka && \
    sed -i 's/tensorflow/##tensorflow-gpu/' requirements.txt && \
    sed -i 's/ --prompt/ --system-site-packages --prompt/' Makefile && \ 
    make install

# acts similar to activating python venv before running
# LC_ALL=C for singularity compatibility
ENV PATH="/medaka/venv/bin:$PATH" \
    LC_ALL=C