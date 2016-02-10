# Dockerfile
#
# Version:          1.0
# Software:         FastQC
# Software Version: 0.11.4
# Description:      Generate statistics on raw data.
# Website:          http://nebula.curie.fr
# Tags:             NGS: QC
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/fastqc:0.11.4 .
# Pull Cmd:         docker pull institutcuriengsintegration/fastqc:0.11.4
# Run Cmd:          docker run --rm -it institutcuriengsintegration/fastqc:0.11.4 <options> <files>
#################################################################

# Set the base image to Ubuntu and python 2.7
FROM ubuntu:14.04
FROM python:2.7

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      default-jre \
      samtools && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/fastqc
WORKDIR /usr/bin/fastqc
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/FastQC.tar.gz -P /usr/bin/fastqc
RUN tar -zxvf FastQC.tar.gz
RUN rm -f FastQC.tar.gz
WORKDIR /usr/bin/fastqc

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
