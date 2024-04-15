FROM python:3.10.12

# for easy upgrade later. ARG variables only persist during image build
ARG SNIFFLESVER=2.3.1

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pysam edlib psutil sniffles==${SNIFFLESVER}

CMD [ "sniffles" ]
