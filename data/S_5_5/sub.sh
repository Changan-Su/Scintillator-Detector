#!/bin/bash
#SBATCH -p amd_a8_384
#SBATCH -N 1
#SBATCH -n 128
source /public1/soft/modules/module.sh
module load cmake/3.30.2-zyq gcc/14.2.0
source /public1/home/a8s001349/soft/geant4-11.2.2/install/bin/geant4.sh
export PATH=/public1/home/a8s001349/soft/geant4-11.2.2/install/bin:$PATH
source /pulbic1/home/a8s001349/soft/g4data.sh
./exampleB1 ../../run4.mac
