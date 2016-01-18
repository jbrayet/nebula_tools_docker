# Dockerfile
#
# Version:          1
# Software:         Macs
# Software Version: 1.4.2
# Description:      This tool allows ChIP-seq peak calling using MACS.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/macs:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/macs:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/macs:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu, python 2.7 and R 3.1.2
FROM ubuntu:14.04
FROM python:2.7
FROM r-base:3.1.2

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/macs
WORKDIR /usr/bin/macs
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/MACS_1.4.2.tar.gz -P /usr/bin/macs
RUN tar -zxvf MACS_1.4.2.tar.gz
RUN rm -f MACS_1.4.2.tar.gz

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
