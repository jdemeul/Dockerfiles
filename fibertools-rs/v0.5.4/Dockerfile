FROM rust:bullseye as builder

# for easy upgrade later. ARG variables only persist during image build
ARG FIBERTOOLSVER=0.5.4

LABEL base.image="rust:bullseye"
LABEL dockerfile.version="1"
LABEL software="fibertools-rs"
LABEL software.version="0.5.4"
LABEL description="NA"
LABEL website="https://github.com/fiberseq/fibertools-rs"
LABEL license="https://github.com/nanoporetech/modkit/blob/master/LICENCE.txt"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y git cmake gcc && \
 apt autoclean && rm -rf /var/lib/apt/lists/*

RUN rustup component add rustfmt
RUN wget https://download.pytorch.org/libtorch/cu118/libtorch-shared-with-deps-2.2.0%2Bcu118.zip && \
    unzip libtorch-shared-with-deps-2.2.0+cu118.zip && rm libtorch-shared-with-deps-2.2.0+cu118.zip && \
    mv libtorch /usr/lib/.

ENV LIBTORCH_CXX11_ABI=0
ENV LIBTORCH=/usr/lib/libtorch
ENV LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
ENV DYLD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH

RUN cargo install --all-features fibertools-rs@${FIBERTOOLSVER}

CMD ["ft"]
# FROM debian:bullseye-slim
# # RUN apt-get update && apt-get install -y extra-runtime-dependencies && rm -rf /var/lib/apt/lists/*
# COPY --from=builder /usr/local/cargo/bin/modkit /usr/local/bin/modkit
# CMD ["modkit"]