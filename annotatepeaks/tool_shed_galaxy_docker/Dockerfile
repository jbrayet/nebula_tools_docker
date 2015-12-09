#Dockerfile
#
# Version:          1
# Software:         AnnotatePeaks
# Software Version: 1.0
# Description:      Tools for annotates peaks with genomic feature (promoter, enhancer, exon, intron, etc.) and creates a .png or .pdf file with distribution.
# Website:          http://nebula.curie.fr
# Tags:             NGS: Peak Annotation
# Provides:         samtools 1.0
# Build Cmd:        docker build --rm -t institutcuriengsintegration/annotatepeaks:1.0 .
# Pull Cmd:         docker pull institutcuriengsintegration/annotatepeaks:1.0
# Run Cmd:          docker run --rm -it institutcuriengsintegration/annotatepeaks:1.0 <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM ubuntu:14.04

RUN apt-get update && \
    apt-get install -y \
      wget && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN 


##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Jocelyn Brayet <jocelyn.brayet@curie.fr>
