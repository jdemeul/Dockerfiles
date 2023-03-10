FROM hkubal/clair3:v1.0.0

ARG RERIOTAG="514e5dd"

LABEL base.image="hkubal/clair3:v1.0.0"
LABEL dockerfile.version="1"
LABEL software="clair3"
LABEL software.version="1.0.0"
LABEL description="Clair3 - Symphonizing pileup and full-alignment for high-performance long-read variant calling with rerio models"
LABEL website="https://github.com/HKU-BAL/Clair3"
LABEL maintainer="Alexandra Pancikova"
LABEL maintainer.email="alexandra.pancikova@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install rerio
RUN git clone https://github.com/nanoporetech/rerio && \
    cd rerio && git checkout ${RERIOTAG} && cd .. && \
    #  Download all models
    python3 rerio/download_model.py --clair3