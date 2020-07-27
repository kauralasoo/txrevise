FROM bioconductor/bioconductor_docker:RELEASE_3_11
LABEL authors="Andreas Vija" \
      description="Docker image containing all requirements for txrevise R package and pipeline"

RUN apt-get update && apt-get install -y libbz2-dev libcurl4 libxml2-dev zlib1g-dev
#RUN apt-get update && apt-get install -y libcurl4-openssl-dev
RUN R -e "install.packages(install.packages('XML', repos = 'http://www.omegahat.net/R'))"
#RUN R -e "install.packages('RCurl')"
RUN R -e "BiocManager::install(c('GenomicFeatures', 'devtools', 'optparse', 'rtracklayer', 'purrrlyr', 'readr', 'tidyr', 'data.table', 'ggplot2'))"
RUN R -e "devtools::install_github('kauralasoo/txrevise')"
