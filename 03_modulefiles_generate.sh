#!/bin/bash
SOFTWARE_PATH=$1
MODULEFILES=$2
mkdir -p $MODULEFILES
cd  $MODULEFILES
for i in `ls $SOFTWARE_PATH`;
do
  if [[ $i == "Leo-xd-as7.8-64-20210910-1" || 
        $i == "Zotero" ||
        $i == "baidunetdisk" ||
        $i == "clash" ||
        $i == "google" ||
        $i == "gromacs" ||
        $i == "icu" ||
        $i == "teamviewer" ||
        $i == "todesk" ||
        $i == "ibm" ||
        $i == "containerd" ||
        $i == "cluster-buster" ||
        $i == "ADFRsuite" ||
        $i == "anaconda2" ||
        $i == "create_cisTarget_databases" ||
        $i == "jupyterhub" ||
        $i == "microsoft" ||
        $i == "NINJA" ||
        $i == "anaconda3" ]]; then 
          continue 
  fi
  mkdir -p $i
  cd $i
  for j in `ls ${SOFTWARE_PATH}/$i`;
  do
    # cat $j && continue && echo "$i/$j already exists!"
    echo "#%Module1.0"                                      > $j
    echo "set  version  $j"                                 >> $j
    echo "set  prefix   ${SOFTWARE_PATH}/$i/$j"               >> $j
    echo "prepend-path  PATH   ${SOFTWARE_PATH}/$i/$j/bin"    >> $j
  done
  
  defaultVer=`ls ${SOFTWARE_PATH}/$i | tail -1`
  echo "#%Module1.0"                              > '.version'
  echo "set ModulesVersion \"$defaultVer\""       >> '.version'

  echo "$i is ready!"

  cd ../
done