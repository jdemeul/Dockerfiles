FROM python:3.9

# for easy upgrade later. ARG variables only persist during image build
ARG SNIFFLESVER=2.0.2

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pysam sniffles==${SNIFFLESVER}

CMD [ "sniffles" ]