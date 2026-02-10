#! /usr/bin/awk -f

{
	OFS = "\t"
}

NR == 1{
	date[$2] = $4
	place[$2] = $3

	print $1,$2
	next
}
{
	if(date[$2] == $4 && place[$2] == $3){
		next
	}else{
		date[$2] = $4
		place[$2] = $3

		print $1,$2
	}
}
