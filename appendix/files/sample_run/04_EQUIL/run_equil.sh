#!/bin/bash
#SBATCH --job-name="WT10_EQUIL"
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --partition=exx96
#SBATCH --mem=10
#SBATCH -B 1:1:1
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL

# cuda
#export CUDA_VISIBLE_DEVICES=`gpu-free | sed 's/,/\n/g' | shuf | head -1`
export CUDA_HOME=/usr/local/n37-cuda-9.2
export PATH=/usr/local/n37-cuda-9.2/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/n37-cuda-9.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/usr/local/n37-cuda-9.2/lib:${LD_LIBRARY_PATH}"

# openmpi
export PATH=/share/apps/CENTOS6/openmpi/1.8.4/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/CENTOS6/openmpi/1.8.4/lib:$LD_LIBRARY_PATH

# amber18
source /share/apps/CENTOS7/amber/amber18/amber.sh


pmemd.cuda -O -i equil0.in -o equil0.out -p ../02_EMIN/WT10.prmtop -c ../03_HEAT/WT10R3_heat.rst -r equil0.rst -ref ../03_HEAT/WT10R3_heat.rst -x mdcrd0 -inf mdinfo0
pmemd.cuda -O -i equil1.in -o equil1.out -p ../02_EMIN/WT10.prmtop -c equil0.rst -r equil1.rst -ref equil0.rst -x mdcrd1 -inf mdinfo1
pmemd.cuda -O -i equil2.in -o equil2.out -p ../02_EMIN/WT10.prmtop -c equil1.rst -r equil2.rst -ref equil1.rst -x mdcrd2 -inf mdinfo2
pmemd.cuda -O -i equil3.in -o equil3.out -p ../02_EMIN/WT10.prmtop -c equil2.rst -r equil3.rst -ref equil2.rst -x mdcrd3 -inf mdinfo3
pmemd.cuda -O -i equil4.in -o equil4.out -p ../02_EMIN/WT10.prmtop -c equil3.rst -r equil4.rst -ref equil3.rst -x mdcrd4 -inf mdinfo4
pmemd.cuda -O -i equil5.in -o equil5.out -p ../02_EMIN/WT10.prmtop -c equil4.rst -r equil5.rst -ref equil4.rst -x mdcrd5 -inf mdinfo5

