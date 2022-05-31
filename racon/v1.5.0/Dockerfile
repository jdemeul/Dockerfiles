# FROM nvidia/cuda:11.4.0-devel-ubuntu20.04
# FROM nvidia/cuda:10.2-devel-ubuntu18.04
FROM ubuntu:focal

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ARG RACONVER=1.5.0

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage python-is-python3 
RUN apt-get update && apt-get install -y build-essential git cmake zlib1g-dev libboost-all-dev doxygen && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install racon
RUN git clone --branch ${RACONVER} --recursive https://github.com/lbcb-sci/racon.git racon

# RUN cd /racon && \
#     mkdir build && \
#     # sed -i 's/racon_params\[-1\]/racon_params\[18\]/' /racon/scripts/racon_wrapper.py && \
#     cmake -DCMAKE_BUILD_TYPE=Release && \
#     # cd build && make || true && \ 
#     # sed -i 's/operator""_format/operator_format/' /racon/vendor/GenomeWorks/3rdparty/spdlog/include/spdlog/fmt/bundled/format.h && \
#     make

ENV PATH="/racon/build/bin/:${PATH}"

CMD [ 'racon' ]