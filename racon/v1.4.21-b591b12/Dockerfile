FROM nvidia/cuda:11.2.1-devel-ubuntu20.04

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG RACONVER=b591b12

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install -y build-essential git cmake zlib1g-dev libboost-all-dev doxygen python-is-python3 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install sniffles binary
RUN git clone --recursive https://github.com/lbcb-sci/racon.git racon && \
    cd racon && \
    mkdir build && \
    cmake -B build -DCMAKE_BUILD_TYPE=Release -Dracon_enable_cuda=ON -Dracon_build_wrapper=ON && \
    cd build && make || true && \ 
    sed -i 's/operator""_format/operator_format/' /racon/vendor/GenomeWorks/3rdparty/spdlog/include/spdlog/fmt/bundled/format.h && \
    make

ENV PATH="/racon/build/bin/:${PATH}"

CMD [ 'racon' ]