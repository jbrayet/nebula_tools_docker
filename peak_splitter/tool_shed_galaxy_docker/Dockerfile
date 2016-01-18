# Dockerfile
#
# Version:          1.0
# Software:         Peak Splitter
# Software Version: 1.0
# Description:      Subdivides peak regions containing more than one site of signal enrichment.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/peaksplitter:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/peaksplitter:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/peaksplitter:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/peaksplitter
WORKDIR /usr/bin/peaksplitter
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/PeakSplitter_v1.tar.gz -P /usr/bin/peaksplitter
RUN tar -zxvf PeakSplitter_v1.tar.gz
RUN rm -f PeakSplitter_v1.tar.gz

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
