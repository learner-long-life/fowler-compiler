#!/bin/sh

mkdir ${1}_${2}
cd ${1}_${2}

cp ../gate .

rm in

echo "1000000" >> in
echo "${2}" >> in
echo "1" >> in
echo "${1}" >> in

#echo "#PBS -m e" >> ${1}_${2}
#echo "#PBS -l walltime=1000:00:00" >> ${1}_${2}
#echo "" >> ${1}_${2}
#echo "cd /home/agf/Solovay/dat4/${1}_${2}" >> ${1}_${2}
#echo "./gate" >> ${1}_${2}

#qsub < ${1}_${2}

bjssub gate
