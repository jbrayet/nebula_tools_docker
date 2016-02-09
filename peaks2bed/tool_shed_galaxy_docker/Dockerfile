# Dockerfile
#
# Version:          1.0
# Software:         Peaks to bed
# Software Version: 1.0
# Description:      This tool creates a .bed file using .peaks output of FindPeaks.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Calling
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/peakstobed:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/peakstobed:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/peakstobed:1.0 <options> <files>
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

RUN mkdir /usr/bin/peakstobed
WORKDIR /usr/bin/peakstobed
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/peaks2bed/tool_shed_galaxy_script/peaks2bed.pl -P /usr/bin/peakstobed
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/create_annotation_files.sh -P /usr/bin
WORKDIR /usr/bin/peakstobed

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
