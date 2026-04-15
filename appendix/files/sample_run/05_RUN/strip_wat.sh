#!/bin/bash
#SBATCH --job-name="ISOFORM_10FULL_TRAJ_STRIP"
#SBATCH --output=FULL_traj_out
#SBATCH --error=FULL_traj_err
#SBATCH --partition=exx96
#SBATCH -B 1:1:1
#SBATCH -N 1
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

cpptraj -i wat.in > strip_wat.log
