# Dockerfile
#
# Version:          1.0
# Software:         AhoPro
# Software Version: 1.3
# Description:      AhoPro is an exact p-value calculator for multiple occurrences of multiples motifs.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Motif Discovery
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/ahopro:1.3 .
# Pull Cmd:         docker pull institutcuriengsintegration/ahopro:1.3
# Run Cmd:          docker run --rm -it institutcuriengsintegration/ahopro:1.3 <options> <files>
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
      libpng12-dev \
      ruby-full \
      imagemagick \
      libmagickwand-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/ahopro
WORKDIR /usr/bin/ahopro
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/AhoPro.1.3.tar.gz -P /usr/bin/ahopro
RUN tar -zxvf AhoPro.1.3.tar.gz
RUN rm -f AhoPro.1.3.tar.gz
WORKDIR /usr/bin/ahopro

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
