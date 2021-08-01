FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG PYCOQCVER=2.5.2

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pycoQC==${PYCOQCVER}

COPY . .

CMD [ "pycoQC" ]