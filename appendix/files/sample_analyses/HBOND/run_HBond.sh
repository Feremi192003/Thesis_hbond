#!/bin/bash
#SBATCH --job-name="I10_WT_HBond"
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --partition=mw256fd
#SBATCH -B 1:1:1
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH --mail-type=BEGIN,END,FAIL

# amber18
source /share/apps/CENTOS7/amber/amber18/amber.sh

cpptraj -i HBond.in > cpptraj_HBond.log
grep -E 'DA_|DC_|DT_|DG_' avg.dat > i10_INTERACT.dat
grep -E 'ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL' i10_INTERACT.dat > WT_i10_DNA_RES_INTERACT.dat
