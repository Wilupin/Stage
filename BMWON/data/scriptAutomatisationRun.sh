#!/bin/bash

echo -e "\nScript bash d'automatisation des calculs BMW pour O(2) en dimension 2!\n"

if [ ! -e data3 ]
then mkdir data3
fi
cd data3

var_rini=0.011
for iii in `seq 0 8`;
do
	var_r=$(echo "$var_rini + $iii*0.001"| bc -l)
	echo "$var_r"
	
	nomdossier="KT$var_r"

	if [ ! -e $nomdossier ]
	then 
		mkdir $nomdossier
		cd $nomdossier
		mkdir data
		echo "d n alpha u0 rmin rmax rht rlt dichoOnOff tmax ompOnOff numberOfThreads">>data.ini
		echo "2.0 2.0 2.0 0.001 -0$var_r -0$var_r -0.4 -1.6 0 300000 1 6">>data.ini
		
		#icpc -fopenmp ../../*.cpp -Ofast -o O2d2
		cp ../KT.01/O2d2 .

		echo -e "#!/bin/bash\n#SBATCH --job-name=O2d2$var_r\n#SBATCH --output=./test.out\n#SBATCH --error=./test.err\n#SBATCH --mincpus=6\n#SBATCH --cpus-per-task=6\n#SBATCH --exclude=ember\nsrun nice -n 19 O2d2">>RunSlurm   
		sbatch RunSlurm

		cd ../
	fi
done
