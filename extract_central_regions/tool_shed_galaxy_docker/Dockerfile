# Dockerfile
#
# Version:          1
# Software:         Extract central regions
# Software Version: 1.0
# Description:      This tool extract central regions of ChIP-seq peaks.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Extract central regions
# Provides:         Nebula
# Build Cmd:        docker build --rm -t institutcuriengsintegration/extractcentralregions:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/extractcentralregions:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/extractcentralregions:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

#############

RUN apt-get update && \
    apt-get install -y \
      build-essential && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /usr/bin/ExtractCentralRegions
RUN wget https://github.com/jbrayet/nebula_tools_docker/raw/master/extract_central_regions/tool_shed_galaxy_script.tar.gz -P /usr/bin/ExtractCentralRegions
WORKDIR /usr/bin/ExtractCentralRegions
RUN tar -xvf tool_shed_galaxy_script.tar.gz && rm tool_shed_galaxy_script.tar.gz
RUN mv tool_shed_galaxy_script/* . && rm -rf tool_shed_galaxy_script extractCentralRegions_wrapper* tool_dependencies.xml

WORKDIR /usr/bin/ExtractCentralRegions

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
