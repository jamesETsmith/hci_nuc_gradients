#!/bin/bash
######################################################################
#                                                                    #
# Author: James E. T. Smith <james.smith9113@gmail.com>              #
# Date: 3/5/2021                                                     #
#                                                                    #
# Sample usage:                                                      #
#   $ ./custom_submit.sh <param> <multiplicity> <species>            #
#                                                                    #
######################################################################

#################
# User settings #
#################
PARAM=$1
MULT=$2
SPECIES=$3

USER_NODES=1
USER_TIME="36:00:00"
SPACK_ENV="rusty-broadwell"
USER_JOB_NAME="${PARAM}_${MULT}_${SPECIES}"
USER_OMP_NUM_THREADS=28
USER_VERBOSE=true
USER_OUTPUT="slurm.out"
USER_TMP="/scratch"
USER_PARTITION="ccq" # Other options: "genx", "ccq", "mem", "bnl", "bnlx", "gpu"


######################
# Print user options #
######################
if ${USER_VERBOSE}
then
    echo "General Settings"
    echo "================"
    echo "Activating Conda Env = ${USER_CONDA_ENV}"
    echo "Wall time limit = ${USER_TIME}"
    echo "Job Name = ${USER_JOB_NAME}"
    echo "Partition = ${USER_PARTITION}"
    echo "Using tmp directory = ${USER_TMP}"
    echo "Number of node(s) = ${USER_NODES}"
    echo "Setting OMP_NUM_THREADS = ${USER_OMP_NUM_THREADS}"
    echo "Slurm output file = ${USER_OUTPUT}"
    echo ""
fi


############################################
# Here is the constructed slurm batch file #
############################################

echo "#!/bin/bash
#SBATCH --output ${USER_OUTPUT}
#SBATCH --job-name ${USER_JOB_NAME}
#SBATCH --time=${USER_TIME}
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --partition=${USER_PARTITION}
#SBATCH --constraint=broadwell

export OMP_NUM_THREADS=${USER_OMP_NUM_THREADS}
export TMP=${USER_TMP}
export TEMP=${USER_TMP}
export TMPDIR=${USER_TMP}

#################
### Node Info ###
#################

echo "Processor Info"
echo "=============="
lscpu
cat /sys/devices/cpu/caps/pmu_name
echo ""

echo "Memory Info"
echo "==========="
head -n 2 /proc/meminfo
echo ""


######################
### Load Spack Env ###
######################
source ~/apps/spack/share/spack/setup-env.sh
spack env activate ${SPACK_ENV} 
source ~/apps/pyscf/pyscf_env.sh
echo ""
echo ""

########################
### Run Calculations ###
########################

../../run_pyscf.sh ${MULT} ${SPECIES}

" > _slurm_batch.sh

sbatch _slurm_batch.sh
rm _slurm_batch.sh # Comment out to see record of slurm job
