#!/bin/bash

echo -e "\nScript bash d'automatisation du script gnuplot!\n"

echo "set title 'Visualisation Flot Eta'">gnuscript2

var_rini=0.01

var_r=$(echo "$var_rini"| bc -l)
	
nomdossier="KT$var_r"
if [ $nomdossier ]
then 
	echo "$nomdossier"
	echo -e "\nplot '$nomdossier/data/w0eta' u 1:4 w l, \ ">>gnuscript2 
fi


for iii in `seq 1 9`;
do
	var_r=$(echo "$var_rini + $iii*0.001"| bc -l)
	
	nomdossier="KT$var_r"

	if [ $nomdossier ]
	then 
		echo "$nomdossier"
		echo -e "'$nomdossier/data/w0eta' u 1:4 w l, \ ">>gnuscript2 
	fi
done

var_r=$(echo "$var_rini+  10*0.001"| bc -l)
	
nomdossier="KT$var_r"
if [ $nomdossier ]
then 
	echo "$nomdossier"
	echo -e "'$nomdossier/data/w0eta' u 1:4 w l \n">>gnuscript2
fi

echo -e "pause -1">>gnuscript2
