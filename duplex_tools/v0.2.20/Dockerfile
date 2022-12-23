FROM python:3.10

# for easy upgrade later. ARG variables only persist during image build
ARG DUPLEXTOOLSVER=0.2.20

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir duplex_tools==${DUPLEXTOOLSVER}
