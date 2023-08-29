FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG SINTOVER=0.9.0

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir sinto==${SINTOVER}

CMD [ "sinto" ]