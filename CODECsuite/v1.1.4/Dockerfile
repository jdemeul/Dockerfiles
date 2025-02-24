FROM ubuntu:18.04 AS build

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update

RUN apt-get install -y autoconf --no-install-recommends \
    git \
    gcc-7=7.3.0-16ubuntu3 \
    gcc-7-base=7.3.0-16ubuntu3 \
    libstdc++-7-dev=7.3.0-16ubuntu3 \
    g++-7=7.3.0-16ubuntu3 \
    cpp-7=7.3.0-16ubuntu3 \
    libgcc-7-dev=7.3.0-16ubuntu3 \
    libasan4=7.3.0-16ubuntu3 \
    libubsan0=7.3.0-16ubuntu3 \
    libcilkrts5=7.3.0-16ubuntu3 \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libssl-dev \
    wget \
    ca-certificates

RUN apt-get install -y --no-install-recommends build-essential && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*
    
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.3/cmake-3.18.3.tar.gz && \
    tar -xzf cmake-3.18.3.tar.gz && \
    cd cmake-3.18.3 && \
    ./bootstrap && \
    make && \
    make install && \
    rm -rf /cmake-3.18.3 cmake-3.18.3.tar.gz

RUN git clone --recursive https://github.com/broadinstitute/CODECsuite.git
RUN cd CODECsuite && mkdir build && cd build && cmake .. && make
RUN mv /CODECsuite/build/codec /usr/bin/codec && rm -rf /CODECsuite

CMD [ "codec" ]

# FROM scratch
# COPY --from=build /usr/bin/codec /
# ENTRYPOINT ["/codec"]