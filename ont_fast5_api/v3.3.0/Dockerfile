FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG ONTFAST5APIVER=3.3.0

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir ont-fast5-api==${ONTFAST5APIVER}