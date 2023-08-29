FROM rust:bullseye as builder

# for easy upgrade later. ARG variables only persist during image build
ARG MODKITVER=0.1.5

LABEL base.image="rust:bullseye"
LABEL dockerfile.version="1"
LABEL software="modkit"
LABEL software.version="0.1.5"
LABEL description="NA"
LABEL website="https://github.com/nanoporetech/modkit/"
LABEL license="https://github.com/nanoporetech/modkit/blob/master/LICENCE.txt"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y git cmake && \
 apt autoclean && rm -rf /var/lib/apt/lists/*

# install samtools, make /data
RUN git clone --depth 1 --branch v${MODKITVER} https://github.com/nanoporetech/modkit.git && \
    cd modkit && cargo install --config net.git-fetch-with-cli=true --path .

FROM debian:bullseye-slim
# RUN apt-get update && apt-get install -y extra-runtime-dependencies && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/cargo/bin/modkit /usr/local/bin/modkit
CMD ["modkit"]