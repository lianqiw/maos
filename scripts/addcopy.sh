files=$(ls ../maos/*.[ch] ../sys/*.[ch] ../lib/*.[ch] ../math/*.[ch] ../skyc/*.[ch] ../tools/*.[ch] ../mex/*.[ch] ../cuda/*.[ch] ../cuda/*.cu)

for file in $files;do
    if grep -q Copyright $file;then
	line=$(grep -m 1 -n www.gnu.org $file |cut -d ':' -f 1)
	echo line=$line
	line=$((line+1)) 
	if [ $line -ne 17 ];then
	    echo skip $file: wrong format
	else
	    sed -i "1,${line}d" $file || continue
	fi
    fi

    if ! grep -q Copyright $file; then
	echo Updated copyright information from $file
	cat copy >>$file.new 
	echo     >>$file.new
	cat $file>>$file.new
	mv $file.new $file
    fi
done

