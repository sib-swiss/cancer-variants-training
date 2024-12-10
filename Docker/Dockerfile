FROM linuxserver/code-server:4.95.3

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget bc tree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install mamba
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
RUN /opt/miniforge3/condabin/mamba init

# Put conda in path so we can use conda activate
ENV PATH=/opt/miniforge3/condabin/:$PATH


# intall required software with conda
COPY environment.yml ./

RUN mamba env create \
    -n ngs-tools \
    -f environment.yml

RUN . ~/.bashrc && conda activate ngs-tools && vep_install -a p --PLUGINS all

COPY settings.json /config/data/User/
COPY keybindings.json /config/data/User/
RUN mkdir -p /config/project
