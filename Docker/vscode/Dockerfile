FROM linuxserver/code-server:4.106.2

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget bc tree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install mamba
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
RUN rm Miniforge3-$(uname)-$(uname -m).sh

# Put conda in path so we can use conda activate
ENV PATH="/opt/miniforge3/bin:$PATH"
RUN mamba --version

# intall required software with conda
COPY environment.yml ./

RUN mamba env create \
    -f environment.yml -y && mamba clean -a -y

# activate the shell
SHELL ["/bin/bash", "-c"]

# Activate environment and install additional packages
RUN source /opt/miniforge3/etc/profile.d/conda.sh && \
    eval "$(mamba shell hook --shell bash)" && \
    mamba activate ngs-tools && \
    mamba remove perl-bio-samtools samtools ensembl-vep --force -y || true && \
    mamba install -c bioconda -c conda-forge ensembl-vep=115.0 samtools=1.22.1 -y && \
    mamba clean -t -y

# activate in the env
RUN source /opt/miniforge3/etc/profile.d/conda.sh && \
    eval "$(mamba shell hook --shell bash)" && \
    mamba activate ngs-tools && \
    vep_install -a p --PLUGINS all

# Settings
COPY settings.json /config/data/User/
COPY keybindings.json /config/data/User/

# Make the folder for the course data
RUN mkdir -p /config/project

# Set default command to activate environment
RUN echo "source /opt/miniforge3/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "conda activate ngs-tools" >> /etc/bash.bashrc && \
    echo "source /opt/miniforge3/etc/profile.d/conda.sh" >> /root/.bashrc && \
    echo "conda activate ngs-tools" >> /root/.bashrc