#!/bin/bash
cd ${HOME}/Logs
SZ=$(du -s /Applications/Tec100/tlm/tlm.log | cut -f 1)
SZBAK=$(du -s tlm.log.bak | cut -f 1)
if [ $SZ -lt $SZBAK ]; then
  let i=0
  while [ -e tlm.log.${i}.gz ]; do
     let i=$i+1
  done
  echo $i
  mv tlm.log.bak tlm.log.${i}
  gzip tlm.log.${i}
fi
cp /Applications/Tec100/tlm/tlm.log tlm.log.bak
