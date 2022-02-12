FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build time
ARG BONITOVER="0.5.1"

LABEL base.image="ubuntu:focal"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y git python3-pip python3-venv && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install taiyaki
RUN git clone --branch v0.5.1 https://github.com/nanoporetech/bonito.git && \
    cd bonito && cd ..

ENV VIRTUAL_ENV="/bonito/venv3"
RUN cd /bonito && python3 -m venv venv3
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN export CUDA_VERSION="113" && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -f https://download.pytorch.org/whl/torch_stable.html -e bonito/

CMD [ "bonito" ]