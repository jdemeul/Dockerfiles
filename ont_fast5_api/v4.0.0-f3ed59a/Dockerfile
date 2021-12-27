FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG ONTFAST5APIVER=f3ed59a

WORKDIR /usr/src/app

RUN git clone https://github.com/nanoporetech/ont_fast5_api && \
    cd ont_fast5_api && git checkout ${ONTFAST5APIVER} && \
    pip install --no-cache-dir --upgrade pip && \
    pip install .