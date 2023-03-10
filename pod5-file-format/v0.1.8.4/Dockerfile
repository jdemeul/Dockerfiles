FROM python:3.10

# for easy upgrade later. ARG variables only persist during image build
ARG POD5TOOLSVER=0.1.8

WORKDIR /usr/src/app

RUN wget https://github.com/nanoporetech/pod5-file-format/releases/download/${POD5TOOLSVER}/lib_pod5-${POD5TOOLSVER}-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl && pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir lib_pod5-${POD5TOOLSVER}-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl && rm lib_pod5-${POD5TOOLSVER}-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
