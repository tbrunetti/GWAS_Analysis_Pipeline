#!/bin/bash

#SBATCH --time=2880
#SBATCH --ntasks=5
#SBATCH --mem=80000
#SBATCH --job-name=prototyping_pipeline
#SBATCH --output=prototyping_pipeline.log
#SBATCH -p longrun

eval "#SBATCH --account=TICR=${USER}"

module load <python location>
module load <R location>
module load <torque location>
module load <libpng location>
module load <freetype location>
module load < gcc location>

cd <into your virtual environment>
source bin/activate
cd <into cloned repository>



chunky run run_GWAS_analysis_pipeline.py <insert all parameters and arguments here>
