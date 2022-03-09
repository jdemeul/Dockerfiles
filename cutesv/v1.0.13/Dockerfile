FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG CUTESVVER=1.0.13

WORKDIR /usr/src/app

# needed older version of setuptools installed to avoid error in PyVCF setup command: use_2to3 is invalid.
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir setuptools==57.5.0 && \
    pip install --no-cache-dir cuteSV==${CUTESVVER}

CMD [ "cuteSV" ]