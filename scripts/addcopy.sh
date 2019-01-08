#files=$(ls ../maos/*.[ch] ../sys/*.[ch] ../lib/*.[ch] ../math/*.[ch] ../skyc/*.[ch] ../tools/*.[ch] ../mex/*.[ch] ../cuda/*.[ch] ../cuda/*.cu)
files=$(ls ../*/*.[ch] ../cuda/*.cu ../*/*.[ch]pp)
	
for file in $files;do
    if grep -q Copyright $file;then
	line=$(grep -m 1 -n www.gnu.org $file |cut -d ':' -f 1)
	echo line=$line
	line=$((line+1)) 
	if [ $line -ne 17 ];then
	    echo skip $file: wrong format
	else
	    if sed -i "1,${line}d" $file ;then
		echo Deleted copyright information from $file
	    fi
	fi
    fi

    if ! grep -q lianqiw-at-tmt-dot-org $file; then
	echo Updated copyright information from $file
	cat copy >$file.new 
	echo     >>$file.new
	cat $file>>$file.new
	mv $file.new $file
    fi
done

