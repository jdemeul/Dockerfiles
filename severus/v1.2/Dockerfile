FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]
ARG SEVERUSVER=1.2

RUN apt-get update && \
    apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip && \
    apt-get install --no-install-recommends -y lsb-release && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN cd /opt && \
    git clone --branch ${SEVERUSVER} https://github.com/KolmogorovLab/Severus

RUN cd /opt/Severus && \
    conda env create --name severus_env --file environment.yml
    # conda activate severus_env
    # ./severus.py 

RUN echo "source activate severus_env" > ~/.bashrc
ENV PATH=/opt/conda/envs/severus_env/bin:$PATH
# ENV PATH=/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH


# RUN conda install -y python=3.6 && \
#     conda install -y -c conda-forge -c bioconda numpy scipy pyvcf pystan==2.17.1.0 pysam pyfaidx

# RUN conda install severus -y

# RUN cd /opt && \
#     git clone https://github.com/ruanjue/bsalign.git && \
#     cd bsalign && make

# ENV PATH=/opt/bsalgin/:$PATH

# git clone https://github.com/KolmogorovLab/Severus
# cd Severus
# conda env create --name severus_env --file environment.yml
# conda activate severus_env
# ./severus.py