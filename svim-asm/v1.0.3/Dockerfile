FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG SVIMASMVER=1.0.3

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir svim-asm==${SVIMASMVER}

CMD [ "svim-asm" ]