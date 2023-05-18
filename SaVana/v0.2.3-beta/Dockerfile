FROM python:3.9.6

# for easy upgrade later. ARG variables only persist during image build
ARG SAVANAVER=0.2.3

RUN apt-get update && apt-get install --no-install-recommends -y \
    git && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/app

RUN git clone -b 0.2.3 https://github.com/cortes-ciriano-lab/savana.git && \
    cd savana && pip install -r requirements.txt && \
    python3 -m pip install . -vv

