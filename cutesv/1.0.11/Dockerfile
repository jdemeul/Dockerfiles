FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG CUTESVVER=1.0.11

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir cuteSV==${CUTESVVER}

CMD [ "cuteSV" ]