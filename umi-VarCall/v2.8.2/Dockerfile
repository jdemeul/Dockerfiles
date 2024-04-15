FROM ubuntu:focal

# for easy upgrade later. ARG variables only persist during image build
ARG WHATSHAPVER=1.1

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="WhatsHap"
LABEL software.version="1.1"
LABEL description="lra is a sequence alignment program that aligns long reads."
LABEL website="https://whatshap.readthedocs.io/en/latest/index.html"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 build-essential \
 python3-dev \
 python3-pip \
 git \
 bwa && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir cython msgpack==0.5.6 scipy pyfaidx pysam statsmodels pyfasta

RUN cd /opt && \
    git clone --recursive https://gitlab.com/vincent-sater/umi-varcal

RUN cd /opt/umi-varcal/functions && \
    rm -f functions.so functions.c functions.cpython-35m-x86_64-linux-gnu.so && \
    python3 setup.py build_ext --inplace

# ENTRYPOINT [ "python3 /opt/umi-varcal/umi-varcal.py" ]

# RUN conda create -y -n {env} python=3.3
# RUN echo "source activate {env}" >/.bashrc
# ENV PATH /opt/conda/envs/{env}/bin:$PATH

# RUN conda install -c bioconda cython msgpack==0.5.6 scipy -y
# RUN conda install -c bioconda pyfaidx pysam -y

# RUN cd /opt && \
#     git clone --recursive https://gitlab.com/vincent-sater/umi-varcal

# RUN cd /opt/umi-varcall/functions && \
#     rm functions.so functions.c functions.cpython-35m-x86_64-linux-gnu.so && \
#     python3 setup.py build_ext --inplace

# # ENV PATH=/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH
# # ENV PATH=/opt/FusionSeeker/:$PATH

# # RUN conda install -y python=3.6 && \
# #     conda install -y -c conda-forge -c bioconda numpy scipy pyvcf pystan==2.17.1.0 pysam pyfaidx

# # RUN cd /opt && \
# #     git clone https://github.com/ruanjue/bsalign.git && \
# #     cd bsalign && make

# ENV PATH=/opt/umi-varcal/:$PATH

