FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget git && \
    # apt-get install --no-install-recommends -y lsb-release && \
    apt-get -qq -y autoremove && \
    apt-get clean

# Install your environment.yaml deps into base env 
# Uncomment once you are ready to start productionizing the image 
# COPY environment.yaml /tmp 
# RUN . /opt/conda/bin/activate && conda env update --name base --file /tmp/environment.yaml 
RUN . /opt/conda/bin/activate && \
    conda install -y python=3.10.14 && \
    conda install -y --channel conda-forge --channel bioconda ncurses bedtools htslib samtools minimap2=2.24

RUN . /opt/conda/bin/activate && \
    conda install -y -c conda-forge -c bioconda pybedtools

RUN . /opt/conda/bin/activate && \
    conda install -y -c conda-forge -c bioconda flair==2.0.0
# conda activate flair

RUN cd /opt && \
    git clone https://github.com/cafelton/FLAIR-fusion-v2.git


# Install your softwares 
# Uncomment once you have software to install 
# RUN mkdir /app 
WORKDIR /opt/FLAIR-fusion-v2/
# COPY * ./ 

# Run the software in conda base environment 
# CMD ["/opt/miniconda/bin/conda", "run", "--no-capture-output", "-n", "base", "python" "app.py"]

# ENV PATH=/opt/vartrix:/opt:$PATH
