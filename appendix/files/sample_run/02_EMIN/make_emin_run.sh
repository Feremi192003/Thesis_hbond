cat << EOF > run_emin.sh
#!/bin/bash
#SBATCH --job-name="WT10R1_EMIN"
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --partition=exx96
#SBATCH -B 1:1:1
#SBATCH -N 1

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

echo "/home/apps/CENTOS7/amber/amber18/bin/sander -O -i emin1.in -p WT10.prmtop -c WT10.inpcrd -r WT10R3.rst -ref WT10.inpcrd -o WT10_emin1.out" >> run_emin.sh 
echo "/home/apps/CENTOS7/amber/amber18/bin/sander -O -i emin2.in -p WT10.prmtop -c WT10R3.rst -r WT10R3_emin2.rst -ref WT10R3.rst -o WT10R3_emin2.out " >> run_emin.sh

for i in $(seq 3 1 12); do
	echo " /home/apps/CENTOS7/amber/amber18/bin/sander -O -i emin${i}.in -p WT10.prmtop -c WT10R3_emin$((i-1)).rst -r WT10R3_emin${i}.rst -ref WT10R3_emin$((i-1)).rst -o WT10_emin_${i}.out" >> run_emin.sh

done
