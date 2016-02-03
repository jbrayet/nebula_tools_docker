# Dockerfile
#
# Version:          1
# Software:         ChIP QC
# Software Version: 1.0
# Description:      This tool creates a .png (or .pdf) file with distribution of peaks heights for ChIP and Control sample.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/chipqc:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/chipqc:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/chipqc:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and R 3.1.2
FROM ubuntu:14.04
FROM r-base:3.1.2

#############

RUN apt-get update && \
    apt-get install -y \
      wget \
      build-essential \
      gfortran \
      libxt-dev \
      xorg-dev \
      libpango1.0-dev \
      libcairo-dev \
      ca-certificates && \   
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/ChIPQC
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/chipqc/tool_shed_galaxy_script.tar.gz -P /usr/bin/ChIPQC
WORKDIR /usr/bin/ChIPQC
RUN tar -xvf tool_shed_galaxy_script.tar.gz && rm tool_shed_galaxy_script.tar.gz
RUN mv tool_shed_galaxy_script/* . && rm -rf tool_shed_galaxy_script annotatePeaks_wrapper* tool_dependencies.xml

WORKDIR /usr/bin/ChIPQC

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
