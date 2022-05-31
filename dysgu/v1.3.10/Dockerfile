FROM python:3.9

USER root

ARG DYSGUVER=1.3.10

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --upgrade pip && \
    pip install numpy && \
    pip install dysgu==${DYSGUVER}

CMD ["/bin/sh"]