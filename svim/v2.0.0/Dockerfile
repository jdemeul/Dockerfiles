FROM python:3

# for easy upgrade later. ARG variables only persist during image build
ARG SVIMVER=2.0.0

WORKDIR /usr/src/app

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir svim==${SVIMVER}

ENV MPLCONFIGDIR=/var/cache/matplotlib

CMD [ "svim" ]