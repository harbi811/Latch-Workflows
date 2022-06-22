FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/wf-base:fbe8-main

# Or use managed library distributions through the container OS's package
# manager.

RUN apt-get update && \
    apt-get install -y build-essential  && \
    apt-get install -y wget &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path 
ENV PATH=$CONDA_DIR/bin:$PATH

# install bwa, samtools,picard, gatk
RUN conda install -c bioconda -y bwa &&\
    conda install -c bioconda -y picard &&\
    conda install -c bioconda -y gatk4 

RUN LANG=en_US.utf-8 &&\
    LC_ALL=en_US.utf-8

RUN apt-get update -y &&\
    apt-get install -y autoconf samtools

# You can use local data to construct your workflow image.  
COPY data /root

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
RUN  sed -i 's/latch/wf/g' flytekit.config
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
