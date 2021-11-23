while read p; do
  echo "$p"
  cp INCAR $p;
  cd $p;
  qsub /u/project/ESS/lstixrud/jd848/metad/pvh/inputs/sub_vasp.sh;
  cd ..;
done <tmp