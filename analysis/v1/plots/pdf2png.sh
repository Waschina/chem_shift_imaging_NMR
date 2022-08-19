pdfs=($(ls | grep ".pdf"))

for i in "${pdfs[@]}"
do
	outp="${i/.pdf/}"
	echo $i
	echo $outp
	
	pdftoppm $i $outp -png -rx 300 -ry 300
done

