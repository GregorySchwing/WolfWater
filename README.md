# WolfWater
Wolf parameter estimation along the vapor-liquid coexistence curve

Running this on the WSU grid environment requires pulling singularity containers before running the workflow.
Make a directory, accessible from compute nodes, for storing the images
$ mkdir ~/singularity_cache
Some grid environments, WSU for example, require using the module to build images.
$ module load singularity
Pull the singularity images individually
$ singularity pull  --name go2432-gomc-cpu.img docker://go2432/gomc:cpu
$ singularity pull  --name go2432-namd-latest.img docker://go2432/namd:latest
$ singularity pull  --name go2432-mosdef-gomc-latest.img docker://go2432/mosdef-gomc:latest
$ singularity pull  --name go2432-scikit-optimize-latest.img docker://go2432/scikit-optimize:latest
$ cd ~/WolfWater
The singularity.conf file has to match that of the grid environment in the conda environment.
$ cp /opt/ohpc/pub/libs/singularity/3.4.1/etc/singularity/singularity.conf /wsu/home/go/go24/go2432/mambaforge/envs/nextflow/etc/singularity
A working running script to launch the nextflow pipeline is included.
It sets some environment variables, activates a mamba environment with nextflow and singularity, and runs the pipeline.
$ sbatch headless.sh


