# Drosophila genome assembly paper workflows
Scripts and Dockerfiles for [our paper](https://www.biorxiv.org/) are provided
here. Using these resources, one should be able to exactly reproduce the
assembly workflow and analyses presented in the manuscript. Details on each set
of analyses are found as READMEs in each subfolder.

## Downloading project files

### Data access through NCBI
The reads and genomes produced by this work are available through NCBI from
BioProject PRJNA675888. A supplementary table with accession numbers is provided
in our manuscript. Alternatively, only for species that we sequenced and
uploaded to NCBI, a table can be generated with the [SRA Run
Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA675888). 

### Data access through request
Genome assemblies are still being processed by NCBI and are not yet available
for download. We are temporarily hosting assemblies at [this
link](https://web.stanford.edu/~bkim331/files/genomes/). 

Nanopore raw data (fast5) files are not publicly available but we will provide
them freely upon request. Please email Bernard Kim (email on GitHub profile
page) for access. 

We recommend Rclone, an Rsync like utility for cloud and web storage, for batch
downloads of files. Download it [here](https://rclone.org/downloads/). 

Example for downloading all genome files provided above:
```bash
rclone copy -P --http-url https://web.stanford.edu/~bkim331/files/genomes/ :http: ./
```

Full list of commands: https://rclone.org/commands/

## Setting up containers
We use Docker and Singularity to manage containers but ultimately prefer to
utilize Singularity image files as they are friendlier with the cluster. Docker
requires elevated permissions for certain operations. 

The system for setting up and building images should have Docker and Singularity
(>3.0) installed. We have built images successfully on Linux systems and on
Windows, running Ubuntu 20.04 with the Windows Subsystem for Linux compatibility
layer (WSL2). At the time this was written, Singularity images could not be
built on Mac OS due to compatibility issues, but the images should run. Some
images may not build or work properly without an NVIDIA graphics card. In those
cases you will have to modify the installation to omit the GPU-enabled versions
of certain programs.  

For Nanopore base calling, Racon polishing, and Medaka polishing, an NVIDIA
graphics card of the Pascal (GTX 1000) generation or later should be installed
as well as NVIDIA/CUDA drivers and NVIDIA Docker. 

### Required programs for building images

The following programs must be installed to build images. Once
a Singularity image has been built, it can be transferred to and run
on any system (e.g. the cluster) without elevated privileges.

Docker: https://docs.docker.com/get-docker/  
Singularity: https://sylabs.io/guides/3.0/user-guide/installation.html  
[Overview](https://github.com/NVIDIA/nvidia-docker) for installing Docker with
the NVIDIA Container toolkit with [detailed instructions](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#docker).
NVIDIA drivers* are easily installed on Ubuntu systems with:    
```bash
ubuntu-drivers devices
sudo ubuntu-drivers autoinstall
```
*Optional for GPU-accelerated images.

Also note, when using containers the CUDA libraries do not have to be installed
on the host machine, just the drivers. Use `nvidia-smi` inside a running
container to check that drivers are installed and working and that the graphics
card is detected. 

### Building Docker/Singularity images

Once Docker and Singularity have been installed, an image is built from the
supplied Dockerfiles. Again, note these Dockerfiles contain installation
instructions for libraries/programs of the *exact same versions* used in our
analysis for the sake of reproducibility. Please use current versions of all
programs for better performance.  

From the directory containing the Dockerfile of interest, run the following code
to build the Docker image. For example, to build the Nanopore assembly image:   
```bash
cd ./dockerfiles/assembly

imageName="assembly"
sudo docker build -t ${imageName} .
```  
Once the image is built, a Docker container can be launched with the image. The 
```--gpus all``` argument allows the container to access GPU resources; omit
this if you are not running a GPU image.
```bash
sudo docker run --gpus all -i -t assembly:latest
```
Because the Docker daemon requires root permissions we use Singularity 
to run containers on the cluster. A Singularity image is built from
the Docker image with the following command:
```bash
imageName="assembly"
sudo singularity build ${imageName}.simg docker-daemon://${imageName}:latest
```  
The resulting `assembly.simg` Singularity image is used to run
containers. For example, to start an interactive shell at the current
working directory, use:
```bash
singularity shell --nv assembly.simg
```
Similar to the Docker command above, omit the ```--nv``` argument if GPUs are
not needed for the workflow. Once the Singularity image is built, it can be
copied to a different system and run immediately.
