FROM rocker/shiny:latest

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    git

WORKDIR /srv/shiny-server/

RUN git clone https://github.com/ExpectationsManaged/STRaitRazorOnline.git

RUN R -e "install.packages(c('BiocManager', 'stringi', 'rhandsontable', 'tidyverse', 'shinydashboard', 'tcltk'), dependencies=TRUE)" \
    && R -e "BiocManager::install(c('GenomeInfoDb', 'Biostrings', 'tidyselect'))"

RUN chown -R shiny:shiny /srv/shiny-server

CMD ["/init"]