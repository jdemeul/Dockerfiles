FROM ubuntu:jammy

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 libboost-all-dev \
 libgtest-dev \
 zlib1g-dev \
 gcc \
 g++ \
 wget \
 make \
 googletest \
 ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# ENV GTEST_ROOT=/usr/src/googletest/googletest/
ENV HAS_GTEST=0
# ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/src/googletest/googletest/include/gtest/

# # install samtools, make /data
RUN wget https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.3.0.tar.gz && \
    tar -xzf 1.3.0.tar.gz && \
    rm 1.3.0.tar.gz && \
    cd DRAGMAP-1.3.0 && \
    make && \
    make install

# # set perl locale settings
# ENV LC_ALL=C
