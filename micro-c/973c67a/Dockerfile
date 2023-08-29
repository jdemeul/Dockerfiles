FROM ubuntu:bionic

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get upgrade && \
    apt-get install -y build-essential python3 python3-pip git

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

RUN apt-get -y install bedtools \
    libjpeg-dev \
    libpng-dev \
    zlib1g-dev \
    libcurl4 \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    liblz4-tool \
    curl

RUN python3 -m pip install --upgrade pip

RUN python3 -m pip install \
    pysam \
    tabulate \
    numpy \
    scipy \
    py2bit \
    matplotlib \
    pyBigWig \
    deeptools \
    pandas

RUN python3 -m pip install pairtools

# RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
#     | tar jxf - && \
#     mv bwa-mem2-2.2.1_x64-linux/bwa-mem2 /usr/local/bin && ln -s /usr/local/bin/bwa-mem2 /usr/local/bin/bwa

RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa; make -j $(nproc) && \
    mv bwa qualfa2fq.pl xa2multi.pl /usr/local/bin

RUN curl -L https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
    | tar jxf - && cd samtools-1.11/ && \
    ./configure && make -j $(nproc) && make install && \
    cd htslib-1.11 && make && make install

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"

RUN curl -L https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz \
    | tar xzf - && cd preseq-3*/ && \
    ./configure --enable-hts CPPFLAGS='-I /home/ubuntu/samtools-1.11/htslib-1.11/' LDFLAGS='-L/home/ubuntu/samtools-1.11/htslib-1.11/' && \
    make -j $(nproc) && make install

RUN apt-get autoclean && rm -rf /var/lib/apt/lists/*

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8