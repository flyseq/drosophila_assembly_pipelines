# Drosophila genome assembly paper workflows
Scripts and Dockerfiles for [PAPER] are provided here. Using these resources, one should be able to exactly reproduce the assembly workflow and analyses presented in the manuscript.

## Downloading project files

### Data access through NCBI
The reads and genomes produced by this work are available through NCBI from BioProject PRJNA675888. A supplementary table with accession numbers is provided in our manuscript. Alternatively, only for species that we sequenced and uploaded to NCBI, a table can be generated with the SRA Run Selector:
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA675888

### Data access through request
Genome assemblies are still being processed by NCBI and are not yet available for download. We are temporarily hosting assemblies at the following link:
https://web.stanford.edu/~bkim331/files/genomes/

Nanopore raw data (fast5) files are not publicly available but we will provide them freely upon request. Please email Bernard Kim (email on GitHub profile page) for access.

The easiest way to access project files is to use Rclone, a command line tool for managing cloud storage. Download it here:
https://rclone.org/downloads/

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
The host system must have Docker and Singularity (>3.0) set up for the container images to build properly. Some images may not build or work properly without an NVIDIA graphics card. In those cases you will have to modify the installation to omit the GPU-enabled versions of certain programs. We have built images successfully on Linux systems and on Windows, running the Windows Subsystem for Linux compatibility layer (WSL2). At the time this was written, Singularity images could not be built on Mac OS due to compatibility issues, but the images should run. We have built images successfully on Ubuntu booting natively or installed alongside Windows with the WSL compatibility layer.

For Nanopore base calling, Racon polishing, and Medaka polishing, an NVIDIA
graphics card of the Pascal (GTX 1000) generation or later should be installed
as well as NVIDIA/CUDA drivers and NVIDIA Docker. In Ubuntu, drivers were
installed by running the commands:

```bash
ubuntu-drivers devices
ubuntu-drivers autoinstall
```

[Overview](https://github.com/NVIDIA/nvidia-docker) for installing Docker with
the NVIDIA Container toolkit with [detailed instructions]
(https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#docker).