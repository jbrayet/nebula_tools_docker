# Dockerfile
#
# Version:          1.0
# Software:         ControlSubSet
# Software Version: 1.0
# Description:      This tool does preprocessing of ChIP-seq aligned reads before calling FindPeaks.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/controlsubset:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/controlsubset:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/controlsubset:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/controlsubset
WORKDIR /usr/bin/controlsubset
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/controlsubset/tool_shed_galaxy_docker/createControlSubset.pl -P /usr/bin/controlsubset

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
