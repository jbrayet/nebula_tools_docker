# Dockerfile
#
# Version:          1.0
# Software:         Chipmunk
# Software Version: 6.0
# Description:      De novo motif finding.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Motif Discovery
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/chipmunk:6.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/chipmunk:6.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/chipmunk:6.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

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

RUN mkdir /usr/bin/chipmunk
WORKDIR /usr/bin/chipmunk
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/ChIPMunk_6.0a.tar.gz -P /usr/bin/chipmunk
RUN tar -zxvf ChIPMunk_6.0a.tar.gz
RUN rm -f ChIPMunk_6.0a.tar.gz
WORKDIR /usr/bin/chipmunk

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
