FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG BONITOVER="0.5.3"

LABEL base.image="ubuntu:focal"
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
    pip install --no-cache-dir -f https://download.pytorch.org/whl/torch_stable.html ont-bonito-cuda113==${BONITOVER} protobuf==3.20.1

# install models
RUN bonito download --models

CMD [ "bonito" ]