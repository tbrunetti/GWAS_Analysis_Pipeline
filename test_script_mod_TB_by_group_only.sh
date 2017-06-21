#!/bin/bash

#SBATCH -p defq # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node
#SBATCH -t 0-30:00         # Running time of 1 hours
#SBATCH --share 

eval "#SBATCH --account=TICR=${USER}"

module load R-3.3.1-gcc-6.1.0-tpkh25pakjkfxwqth7kp3gxl4grlf6c2
module load torque
module load gcc/5.1.0

  #split into chunks with 100000 lines
bim=$pathPLINKprefix".bim"
echo $bim
split -d -l 100000 $bim $pathPLINKprefix"_split"
  #reattach the header to each and clean up
for j in $pathPLINKprefix"_split"*;
  do
    echo $j
    out=$j
    echo $out
    /homelink/brunettt/TOOLS/plink --bfile ${pathPLINKprefix} --extract $j --make-bed --out $out
    bed_file=${out}".bed"
    bim_file=${out}".bim"
    fam_file=${out}".fam"
    gds_file=${out}".gds"
    #subset_val=$(echo $j|cut -d "." -f 2)
    Rscript GENESIS_analysis_mod_TB.R $out $pathPLINKprefix $rlib $pcs
 done
