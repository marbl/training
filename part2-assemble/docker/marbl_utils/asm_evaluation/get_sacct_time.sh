#/usr/bin/bash
#calculate time in hours based on verkko log, works with HiC

outf=$1.sacct
rm -f $outf
sacct -j $(grep jobid $1 | grep -oP "(?<=external jobid ')[0-9]+" | tr '\n' ',') --format="JobID,User,CPUTime,CPUTimeRAW,TotalCPU,State,MaxRSS,MaxVMSize,ReqMem" > $outf

awk '{print $5}' $outf |awk -F '[:-]' '{if (NF == 4) {sum += ($1*24 + $2)*3600 + $3*60 + $4} else if (NF==3) {sum += $1*3600 + $2*60 + $3} else if (NF==2) {sum += $1*60 + $2}} END { printf "%.2f\n", sum/3600 }'
