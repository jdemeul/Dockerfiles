FROM ubuntu:jammy

# for easy upgrade later. ARG variables only persist during image build time
ARG SMRTLINKVER="11.0.0.146107"

# metadata
LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="SMRT Link"
LABEL maintainer.email="jonas.demeulemeester@kuleuven.be"

ENV TZ=Europe/Brussels
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install deps and cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y curl git ca-certificates unzip rsync locales && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Set the locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# install SMRT Link
RUN curl -O -L https://downloads.pacbcloud.com/public/software/installers/smrtlink_${SMRTLINKVER}.zip && \
    unzip smrtlink_${SMRTLINKVER}.zip && rm smrtlink_${SMRTLINKVER}.zip

ENV SMRT_ROOT=/opt/pacbio/smrtlink
ENV SMRT_USER=smrtanalysis
RUN chmod 777 /opt

RUN useradd -ms /bin/bash smrtanalysis
USER smrtanalysis

RUN ./smrtlink_${SMRTLINKVER}.run --rootdir $SMRT_ROOT --skip-system-check --skip-userquery --batch --skip-dnscheck
USER root
RUN rm smrtlink_${SMRTLINKVER}.run
USER smrtanalysis

ENV PATH=/opt/pacbio/smrtlink/install/smrtlink-release_11.0.0.146107/smrtcmds/bin/:$PATH