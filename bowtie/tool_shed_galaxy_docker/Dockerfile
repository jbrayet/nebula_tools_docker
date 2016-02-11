# Dockerfile
#
# Version:          1
# Software:         Bowtie
# Software Version: 1.1.2
# Description:      Bowtie_ is a short read aligner designed to be ultrafast and memory-efficient.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Mapping
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/bowtie:1.1.2 .
# Pull Cmd:         docker pull institutcuriengsintegration/bowtie:1.1.2
# Run Cmd:          docker run --rm -it institutcuriengsintegration/bowtie:1.1.2 <options> <files>
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

RUN mkdir /usr/bin/bowtie
WORKDIR /usr/bin/bowtie
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/bowtie-1.1.2.tar.gz -P /usr/bin/bowtie
RUN tar -zxvf bowtie-1.1.2.tar.gz
RUN rm -f bowtie-1.1.2.tar.gz
WORKDIR /usr/bin/bowtie

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
