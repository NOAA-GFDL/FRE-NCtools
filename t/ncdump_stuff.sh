folders=`ls`
for item in $folders[*]}
do

	if [[ $item != *"input"* ]] ; then
  		continue
	fi

	cd $item
	files=`ls`

	for file in $files[*]}
	do
		if [[ -d $file ]]; then
			cd $file
			subdir=`ls`
			for subfile in $subdir[*]}
			do
				if [[ $subfile != *".nc" ]] ; then
  				continue
				fi
				ncdump $subfile > ${subfile}l
				rm $subfile
			done
			cd ..
			continue
		fi

		if [[ $file != *".nc" ]] ; then
  			continue
		fi
		ncdump $file > ${file}l
		rm $file
	done

	cd ..

done

