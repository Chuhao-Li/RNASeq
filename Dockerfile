FROM r-base:4.2.2

MAINTAINER lichuhao<yc27658@um.edu.mo>

# install dependency, libgsl-dev is for gfold. 
RUN apt update && apt install -y fastp bowtie2 subread samtools libgsl-dev python3 pip 

RUN for p in ggplot2 patchwork reshape2 ggdendro argparse BiocManager; do Rscript -e "install.packages('${p}')"; done

# dynamic lib for DESeq2 and sva
RUN apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev xml2 unixodbc-dev

RUN Rscript -e "BiocManager::install('DESeq2'); BiocManager::install('sva');"

RUN pip install jinja2 pandas

COPY --chown=1 gfold /opt/gfold

WORKDIR "/opt/gfold"

RUN make && mv gfold /usr/bin

COPY --chown=1 rnaseq /opt/rnaseq

CMD "/opt/rnaseq/rnaseq.py"

