#!/bin/bash

#SBATCH --time=2880
#SBATCH --ntasks=5
#SBATCH --mem=80000
#SBATCH --job-name=prototyping_pipeline
#SBATCH --output=prototyping_pipeline.log
#SBATCH -p longrun

eval "#SBATCH --account=TICR=${USER}"

cd /gpfs/share/TICR/TICR_Analysis_Pipeline_Prototyping/myVirtENV 
source bin/activate

