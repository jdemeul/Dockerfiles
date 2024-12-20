FROM python:3.10.15

# for easy upgrade later. ARG variables only persist during image build
ARG SNIFFLESVER=2.5.2

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pysam edlib psutil sniffles==${SNIFFLESVER}

CMD [ "sniffles" ]
