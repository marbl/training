FROM ubuntu:latest

COPY commands.sh /commands.sh
COPY marbl_utils /marbl_utils

# Install Miniconda
RUN apt-get update && \
    apt-get install -y wget && \
    apt-get install -y git && \
    apt-get install -y make && \
    apt-get install -y g++ && \
    apt-get install -y zlib1g-dev && \
    apt-get install -y vim && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/conda/bin/conda init && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Set up Conda environment and install packages
RUN /opt/conda/bin/conda create -n myenv python=3.8 && \
    echo "conda activate myenv" >> ~/.bashrc && \
    /opt/conda/bin/conda install -n myenv -c conda-forge -c bioconda -c defaults bedtools mash mashmap minimap2 seqtk winnowmap graphaligner verkko
RUN git clone https://github.com/snurk/gfacpp.git
RUN cd gfacpp && git clone https://github.com/lh3/gfatools.git && make && cd ..
COPY verkko /opt/conda/envs/myenv/bin/
COPY 8-hicPipeline.sm /opt/conda/envs/myenv/lib/verkko/Snakefiles/

# Set the default command to run when the container starts
CMD ["bash"]
