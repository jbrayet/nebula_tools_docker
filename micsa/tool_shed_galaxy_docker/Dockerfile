# Dockerfile
#
# Version:          1.0
# Software:         MICSA
# Software Version: 1.0
# Description:      Motif Identification for Chip-Seq Analysis.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/micsa:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/micsa:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/micsa:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu, java 6 and ruby 2.0
FROM ubuntu:14.04
FROM java:6
FROM ruby:2.0

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential \
      wget \
      samtools && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/micsa
WORKDIR /usr/bin/micsa
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/MICSA.tar.gz -P /usr/bin/micsa
RUN tar -zxvf MICSA.tar.gz
RUN rm -f MICSA.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/ChIPMunk_6.0a.tar.gz -P /usr/bin/micsa
RUN tar -zxvf ChIPMunk_6.0a.tar.gz
RUN rm -f ChIPMunk_6.0a.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_sources/VancouverPackage-4.0.9.2.tar.gz -P /usr/bin/micsa
RUN tar -zxvf VancouverPackage-4.0.9.2.tar.gz
RUN rm -f VancouverPackage-4.0.9.2.tar.gz
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/splitChr.pl -P /usr/bin

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>