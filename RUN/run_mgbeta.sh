#!/bin/bash -l

##SBATCH --account nems
##SBATCH --account da-cpu
#SBATCH --account fv3-cpu
#SBATCH --job-name=run_mgbeta400
##SBATCH --qos batch
#SBATCH -q "debug"
#SBATCH --output=beta400_%j.out
##SBATCH --partition theia
##SBATCH --time=15
##SBATCH --time=8
##SBATCH --time=50
#SBATCH --time=00:05:00
##SBATCH --ntasks=88 --tasks-per-node=4
##SBATCH --ntasks=88  --tasks-per-node=8
##SBATCH --ntasks=88 --tasks-per-node=11
##SBATCH --ntasks=64 --tasks-per-node=16
##SBATCH --ntasks=64 --tasks-per-node=4
##SBATCH --ntasks=384 --tasks-per-node=16
#SBATCH --ntasks=400 --tasks-per-node=16
##SBATCH -D .

set -aeux

module load intel/18.0.5.274
module load impi/2018.4.274


export OMP_NUM_THREADS=1
#export I_MPI_PMI_LIBRARY=/path/to/slurm/pmi/library/libpmi.so
#export I_MPI_PROCESS_MANAGER=mpd

rundir=/scratch1/NCEPDEV/stmp2/Miodrag.Rancic/beta_loc
exedir=/scratch1/NCEPDEV/da/Miodrag.Rancic/EnsLoc/EXE
cov_exe=$exedir/beta_loc.exe

#mkdir $rundir
cd $rundir
rm -f *.dat f
rm -f fort* 
rm -f work* 
rm -f *nml *csv
rm -f *csv
cp /scratch1/NCEPDEV/da/Miodrag.Rancic/EnsLoc/RUN/mgbf.nml $rundir/mgbeta.nml

srun $cov_exe

echo "Done!"
