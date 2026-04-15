#!/bin/bash
#SBATCH --job-name="WT_iso10_RMSD"
#SBATCH --output=traj_out
#SBATCH --error=traj_err
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
cd "/zfshomes/fcaballero/isoformProject/new_wild_type_runs/rep3/iso10/05_PROD/06_POST/RMSD"
cpptraj -i rmsd_avg_struct.in > rmsd_avg_struct.log
