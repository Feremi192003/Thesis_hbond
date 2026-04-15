#!/bin/bash
#SBATCH --job-name="WT10_HEAT"
#SBATCH --output=iso2_out
#SBATCH --error=iso2_err
#SBATCH --partition=exx96
#SBATCH -B 1:1:1
#SBATCH -N 1
#SBATCH --mem=5
#SBATCH --mail-type=BEGIN,END,FAIL

# cuda
export CUDA_HOME=/usr/local/n37-cuda-9.2
export PATH=/usr/local/n37-cuda-9.2/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/n37-cuda-9.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/usr/local/n37-cuda-9.2/lib:${LD_LIBRARY_PATH}"

# openmpi
export PATH=/share/apps/CENTOS6/openmpi/1.8.4/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/CENTOS6/openmpi/1.8.4/lib:$LD_LIBRARY_PATH

# amber18
source /share/apps/CENTOS7/amber/amber18/amber.sh

sander.MPI -O -i heat.in -o heat.out -p ../02_EMIN/WT10.prmtop -c  ../02_EMIN/WT10R3_emin12.rst -r WT10R3_heat.rst -ref ../02_EMIN/WT10R3_emin12.rst
