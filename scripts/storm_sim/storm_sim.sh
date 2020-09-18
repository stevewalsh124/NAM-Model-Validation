#!/bin/bash

# Make sure to change the -A line to your allocation

#SBATCH -N 1
#SBATCH --ntasks-per-node=3
#SBATCH -t 12:00:00
#SBATCH -p normal_q
#SBATCH -A Precipit

# To run this in a loop from 1 to 47, run the following on the command line:
# for i in $(seq 1 47); do sbatch storm_fixnug/storm43fix.sh $i; done

# Add modules
module purge
module load singularity

# Change to the directory from which the job was submitted
#cd $PBS_O_WORKDIR

# Run R


low=${1:-1}
seed=${2:-1}

for i in $(seq $low $low)
do
#echo "Starting batch with seed $i now."
#nohup R CMD BATCH "--args seed=$i reps=$reps" spam_mc.R spam_mc_$i.Rout &

echo "Your requested storm to analyze is $i."

echo "$( date ): Starting storm $i"
singularity exec --bind=/groups:/groups /groups/arcsingularity/ood-rstudio-geospatial_3.6.2.sif R CMD BATCH "--args storms_to_eval=$i seed=$seed" /home/walsh124/NAM-Model-Validation/scripts/sim_REML_Hessian.R hurricane${i}fixnug0.Rout 
echo "$( date ): Finished storm $i" 
done

exit
