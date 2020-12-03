# Drosophila genome assembly paper workflows
Scripts and Dockerfiles for (our paper)[https://www.biorxiv.org/] are provided here. Using these resources, one should be able to exactly reproduce the assembly workflow and analyses presented in the manuscript.

## Downloading project files

### Data access through NCBI
The reads and genomes produced by this work are available through NCBI from BioProject PRJNA675888. A supplementary table with accession numbers is provided in our manuscript. Alternatively, only for species that we sequenced and uploaded to NCBI, a table can be generated with the (SRA Run Selector)[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA675888].

### Data access through request
Genome assemblies are still being processed by NCBI and are not yet available for download. We are temporarily hosting assemblies at (this link)[https://web.stanford.edu/~bkim331/files/genomes/].

Nanopore raw data (fast5) files are not publicly available but we will provide them freely upon request. Please email Bernard Kim (email on GitHub profile page) for access.

The easiest way to access project files is to use Rclone, a command line tool for managing cloud storage. Download it (here)[https://rclone.org/downloads/].

To set it up to work with our Box folder, copy the ```rclone``` binary into a 
directory in your ```$PATH```. Then, configure it by running ```rclone config```
and by choosing the following options. 
1. Add a new remote (```n```) 
1. Name it ```box```, or whatever you'd like
1. Choose option ```6``` for Box storage
1. Don't enter anything for Client Id
1. Don't enter anything for the Client Secret
1. Don't enter anything for the config.json location
1. Sub type is ```1``` or ```"user"```
1. Enter ```n``` and don't edit advanced config
1. If you are working on your own machine, enter ```y```, if you are trying to 
set up Rclone on a remote machine (e.g. the cluster), enter ```n``` 
   - If you're setting up on a local machine, an authorization web page should
     automatically appear. If it does not, go to the address provided by Rclone.
   - If you're setting up on a remote machine, you will need to have a copy of
     ```rclone``` installed on your local machine. At this step you will need to
     run the command ```rclone authorize "box"``` on your local machine, 
     authorize using the web page that shows up, then provide the authorization
     key to the remote machine.

Once ```rclone``` has been set up in this manner, downloading a Box directory is
as easy as running:

```bash
rclone copy -P box:100x100/assemblies/repeat_masked/fastas /path/to/local/dir
```

If you have download some files already and want to skip them:
```bash
rclone sync -P box:100x100/assemblies/repeat_masked/fastas /path/to/local/dir
```

Full list of commands: https://rclone.org/commands/

## Setting up containers
We use Docker and Singularity to manage containers but ultimately prefer to utilize Singularity image files as they are friendlier with the cluster. Docker requires elevated permissions for certain operations.

The system for setting up and building images should have Docker and Singularity (>3.0) installed. We have built images successfully on Linux systems and on Windows, running Ubuntu 20.04 with the Windows Subsystem for Linux compatibility layer (WSL2). At the time this was written, Singularity images could not be built on Mac OS due to compatibility issues, but the images should run. Some images may not build or work properly without an NVIDIA graphics card. In those cases you will have to modify the installation to omit the GPU-enabled versions of certain programs. 

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
the NVIDIA Container toolkit with [detailed instructions]
(https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#docker).
NVIDIA drivers* are easily installed on Ubuntu systems with:    
```bash
ubuntu-drivers devices
sudo ubuntu-drivers autoinstall
```
*Optional for GPU-accelerated images.

Also note, when using containers the CUDA libraries do not have to be installed on the host machine, just the drivers. Use `nvidia-smi` inside a running container to check that drivers are installed and working and that the graphics card is detected.

### Building Docker/Singularity images

Once Docker and Singularity have been installed, an image is built from the supplied Dockerfiles. Again, note these Dockerfiles contain installation instructions for libraries/programs of the *exact same versions* used in our analysis for the sake of reproducibility. Please use current versions of all programs for better performance. 

From the directory containing the Dockerfile of interest, run the following code to build the Docker image. For example, to build the Nanopore assembly image:  
```bash
cd ./dockerfiles/assembly

imageName="assembly"
sudo Docker build -t ${imageName} .
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
sudo singularity build ${imageName}.simg \
    docker-daemon://${imageName}:latest
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