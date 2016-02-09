# Dockerfile
#
# Version:          1.0
# Software:         CCAT intervals to bed
# Software Version: 1.0
# Description:      This tool creates a .bed file using .peaks output of CCAT.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/ccattobed:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/ccattobed:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/ccattobed:1.0 <options> <files>
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

RUN mkdir /usr/bin/ccattobed
WORKDIR /usr/bin/ccattobed
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/ccat2bed/tool_shed_galaxy_script/peaks2bed.pl -P /usr/bin/ccattobed
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
WORKDIR /usr/bin/ccattobed

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
