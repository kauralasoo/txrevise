FROM bioconductor/release_core2
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for txrevise R package and pipeline"

RUN R -e "BiocManager::install(c('GenomicFeatures','devtools', 'optparse', 'rtracklayer','purrrlyr'))"
RUN R -e "devtools::install_github('kauralasoo/txrevise')"
