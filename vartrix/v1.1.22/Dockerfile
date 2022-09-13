FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build time
ARG VARTRIXVER="1.1.22"

# metadata
LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="Vartrix"
LABEL software.version="1be3d4c"
LABEL description="Single-Cell Genotyping Tool."
LABEL website="https://github.com/10XGenomics/vartrix"
LABEL maintainer="Jonas Demeulemeester"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y wget ca-certificates && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install shasta binary
RUN cd /opt && \
    wget https://github.com/10XGenomics/vartrix/releases/download/v${VARTRIXVER}/vartrix_linux && \
    mv vartrix_linux vartrix && \
    chmod 777 vartrix && \
    ln -s /opt/vartrix /usr/local/bin/vartrix

CMD [ 'vartrix' ]