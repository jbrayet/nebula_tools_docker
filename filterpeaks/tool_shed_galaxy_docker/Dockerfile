# Dockerfile
#
# Version:          1.0
# Software:         Find Peaks
# Software Version: 4.0.9.2
# Description:      Peak Finder/Analysis application for the ChIP-Seq.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/findpeaks:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/findpeaks:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/findpeaks:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools \
      default-jre && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/findpeaks
WORKDIR /usr/bin/findpeaks
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/VancouverPackage-4.0.9.2.tar.gz -P /usr/bin/findpeaks
RUN tar -zxvf VancouverPackage-4.0.9.2.tar.gz
RUN rm -f VancouverPackage-4.0.9.2.tar.gz
WORKDIR /usr/bin/findpeaks

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
