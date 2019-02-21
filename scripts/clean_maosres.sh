#!/bin/sh
#Clean maos results folder, leaving only Res_${seed}.bin, Res_${seed}.done, *.conf, *.log

for fd in "$@" ;do
    rm -rfv ${fd}/setup ${fd}/Res[a-zA-Z]*.bin ${fd}/Timing_*.bin ${fd}/psd*.bin
done
    
