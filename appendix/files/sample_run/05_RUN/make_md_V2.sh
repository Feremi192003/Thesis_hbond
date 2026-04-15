cat << EOF > run.sh
#!/bin/bash
#SBATCH --job-name="WT10_MD"
#SBATCH --output=md.out
#SBATCH --error=md.err
#SBATCH --partition=exx96
#SBATCH -B 1:1:1
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL

# cuda
export CUDA_HOME=/usr/local/n37-cuda-9.2
export PATH=/usr/local/n37-cuda-9.2/bin:\$PATH
export LD_LIBRARY_PATH=/usr/local/n37-cuda-9.2/lib64:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/usr/local/n37-cuda-9.2/lib:\${LD_LIBRARY_PATH}"

# openmpi
export PATH=/share/apps/CENTOS6/openmpi/1.8.4/bin:\$PATH
export LD_LIBRARY_PATH=/share/apps/CENTOS6/openmpi/1.8.4/lib:\$LD_LIBRARY_PATH

# amber18
source /share/apps/CENTOS7/amber/amber18/amber.sh

export CUDA_VISIBLE_DEVICES=\`gpu-free | sed 's/,/\n/g' | shuf | head -1\`

EOF
 
echo "pmemd.cuda -O -i md.in -p ../02_EMIN/WT10.prmtop -o output/WT10_md1.out -c ../04_EQUIL/equil5.rst -r output/WT10_md1.rst -x output/WT10_md1.trj" >> run.sh
echo "bash strip_concat.sh strip ../02_EMIN/WT10.prmtop output/WT10_md1.trj output/WT10_md1.nc 1" >> run.sh
echo "rm output/WT10_md1.trj" >> run.sh
for i in $(seq 2 1 1100); do
  echo "pmemd.cuda -O -i md.in -p ../02_EMIN/WT10.prmtop -o output/WT10_md${i}.out -c output/WT10_md$(( $i - 1 )).rst -r output/WT10_md${i}.rst -x output/WT10_md${i}.trj" >> run.sh
  echo "bash strip_concat.sh strip ../02_EMIN/WT10.prmtop output/WT10_md${i}.trj output/WT10_md${i}.nc ${i}" >> run.sh
  echo "rm output/WT10_md${i}.trj" >> run.sh
done

