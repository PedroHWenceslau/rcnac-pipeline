FROM rocker/r-ver:4.3.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

RUN R -e "install.packages(c( \
    'pheatmap', \
    'dplyr', \
    'tidyr', \
    'ggplot2', \
    'gridExtra', \
    'readr', \
    'RColorBrewer', \
    'VennDiagram', \
    'optparse' \
    ), repos='https://cloud.r-project.org')"

WORKDIR /pipeline
