FROM python:3.10

# for easy upgrade later. ARG variables only persist during image build
ARG POD5TOOLSVER=0.1.4

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pod5==${POD5TOOLSVER}
