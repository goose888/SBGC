#!/bin/csh

# Tempelate for creating batch scripts
if ( "$2" == "" ) then  
 
  cat > ${1}.sh << EOF1
#!/bin/sh
#SBATCH -n 1 
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p f,d,g

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.${1} > isam_${1}.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

else

  cat > ${1}_${2}.sh << EOF1
#!/bin/sh
#SBATCH -n 1 
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p f,d,g

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.${1} > isam_${1}.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

endif
