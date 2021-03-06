# Dockerfile
#
# Version:          1
# Software:         HMCan
# Software Version: 1.2
# Description:      HMCan detects histone modifications in cancer samples.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/hmcan:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/hmcan:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/hmcan:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and python 2.7
FROM ubuntu:14.04
FROM python:2.7

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/HMCan
WORKDIR /usr/bin/HMCan
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/HMCanV1.20.tar.gz -P /usr/bin/HMCan
RUN tar -zxvf HMCanV1.20.tar.gz
RUN rm -f HMCanV1.20.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/splitChr.pl -P /usr/bin

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
