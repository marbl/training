SCRIPT_PATH=$BASH_SOURCE
SCRIPT_PATH=`dirname $SCRIPT_PATH`

TOTAL_CPU=0
MAX_MEM=0

asm=$(realpath $1)
for step in `ls $asm/batch-scripts/*.out |awk -F "." '{print $(NF-2)}'|sort |uniq`; do
   TSTEP=0
   TMEM=0
   for name in `ls $asm/batch-scripts/*.$step.*out |awk  -F "/" '{print $NF}'`; do
      id=`echo $name |awk -F "." '{print $1}'`
      res=`sh $SCRIPT_PATH/getTimeAcct.sh $id |grep -v CPU`
      TOTAL_CPU=`echo "$TOTAL_CPU $res" |awk '{print $1+$3}'`
      TSTEP=`echo "$TSTEP $res" |awk '{print $1+$3}'`
      TMEM=`echo "$TMEM $res" |awk '{if ($1 > $4) print $1; else print $4}'`
      WALL=`echo "$WALL $res" |awk '{print $1+$2}'`
      MAX_MEM=`echo "$MAX_MEM $res" |awk '{if ($1 > $4) print $1; else print $4}'`
      echo "$step, $id, $name: $res"
   done
   echo "*****************************************"
   echo "Total cpu for $step: $TSTEP max mem $TMEM" 
done

echo "****************************"
echo "Total CPU TOTAL WALL MAX_MEM"
echo "$TOTAL_CPU $WALL $MAX_MEM"
