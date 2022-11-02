FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG POD5TOOLSVER=0.0.23

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pod5-format-tools==${POD5TOOLSVER}
