FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG DUPLEXTOOLSVER=0.2.14

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir duplex_tools==${DUPLEXTOOLSVER}
