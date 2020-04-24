FROM bioconductor/bioconductor_docker:RELEASE_3_10
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for txrevise R package and pipeline"

RUN R -e "BiocManager::install(c('GenomicFeatures','devtools', 'optparse', 'rtracklayer','purrrlyr','readr', 'tidyr'))"
RUN R -e "devtools::install_github('kauralasoo/txrevise')"
