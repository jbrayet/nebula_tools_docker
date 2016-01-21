# Dockerfile
#
# Version:          1
# Software:         CCAT
# Software Version: 3.0
# Description:      This tool allows ChIP-seq peak/region calling using CCAT.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/ccat:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/ccat:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/ccat:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu and python 2.7
FROM ubuntu:14.04
FROM python:2.7

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/ccat
WORKDIR /usr/bin/ccat
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/CCAT3.0.tar.gz -P /usr/bin/ccat
RUN tar -zxvf CCAT3.0.tar.gz
RUN rm -f CCAT3.0.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/splitChr.pl -P /usr/bin
WORKDIR /usr/bin/ccat

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
