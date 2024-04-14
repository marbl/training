This section of the tutorial focuses on what to do once you have an assembly. All the commands and tools required for the tutorial are packaged in a docker container. If you want to install them locally, the list is:

System applications:
- wget
- git
- make
- g++
- zlib

Bioinformatics tools (available through conda via -c conda-forge -c bioconda -c defaults):
- bedtools
- mash
- mashmap
- minimap2
- seqtk
- winnowmap
- GraphAligner
- verkko
  
Bioinformatics tools (only available through github)
- gfacpp (https://github.com/snurk/gfacpp)

The steps and components required to build the docker image are in the subfolder docker. The commands.sh file lists the commands run as part of the tutorial along with some brief documentation. For more detailed usage information consult the original tool support pages. Some commands are too computationally expensive to run in a docker container and so they are commeted out assuming you'd run them on your own environment.

The commands to download and run the tutorial using the existing image/data are:

    curl -L https://s3.amazonaws.com/genomeark/trainingmaterials/AGBTAG2024/assembly-102.docker.tar.gz -o assembly-102.tar.gz
    curl -L https://s3.amazonaws.com/genomeark/trainingmaterials/AGBTAG2024/tutorial_assembly_data.tar.gz -o tutorial.tar.gz

    # load image into your docker dashboard
    docker load -i assembly-102.tar.gz
    tar xvzf tutorial.tar.gz
    cd tutorial

    # launch docker
    docker run -v "`pwd`:/data"  --platform linux/amd64 -it  assembly-102 /bin/bash

Now just follow the commands.sh steps to try it out.

The assembly folder contains a subset of the tutorial.tar.gz and results of some of the commands in the tutorial (to filter rDNA and satellite repeat from the graph). These are useful to look at the assembly in BandageNG without having to run/download the entire tutorial.
