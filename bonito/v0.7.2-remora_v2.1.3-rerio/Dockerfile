FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build time
ARG BONITOVER="0.7.2"
ARG REMORAVER="v2.1.3"

LABEL base.image="ubuntu:jammy"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y git python3-pip python3-venv && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install bonito
ENV VIRTUAL_ENV="/bonito/venv3"
RUN mkdir /bonito && python3 -m venv venv3
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cu113 ont-bonito

# install models
RUN bonito download --models

# install remora
RUN git clone https://github.com/nanoporetech/remora && \
    cd remora && git checkout ${REMORAVER} && cd ..

RUN pip install --no-cache-dir --upgrade pip && \
    pip install onnxruntime-gpu && \
    pip install --no-cache-dir -e remora/

# install rerio
RUN git clone https://github.com/nanoporetech/rerio && \
    rerio/download_model.py --remora

CMD [ "bonito" ]