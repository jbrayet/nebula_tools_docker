# Dockerfile
#
# Version:          1
# Software:         AnnotatePeaks
# Software Version: 1.0
# Description:      Tools for annotates peaks with genomic feature (promoter, enhancer, exon, intron, etc.) and creates a .png or .pdf file with distribution.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Annotation
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/annotatepeaks:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/annotatepeaks:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/annotatepeaks:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and R 3.1.2
FROM ubuntu:14.04
FROM r-base:3.1.2

#############

RUN apt-get update && \
    apt-get install -y \
      wget \
      curl \
      build-essential \
      gfortran \
      libxt-dev \
      xorg-dev \
      libpango1.0-dev \
      libcairo-dev && \   
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/annotatePeaks
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/annotatepeaks/tool_shed_galaxy_script.tar.gz -P /usr/bin/annotatePeaks
WORKDIR /usr/bin/annotatePeaks
RUN tar -xvf tool_shed_galaxy_script.tar.gz && rm tool_shed_galaxy_script.tar.gz
RUN mv tool_shed_galaxy_script/* . && rm -rf tool_shed_galaxy_script annotatePeaks_wrapper* tool_dependencies.xml

WORKDIR /usr/bin/annotatePeaks

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
