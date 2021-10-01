#Note, remove "" after sed -i in linux
#minus sign must be the last or first inside a []. no need to escape ><
var="[()0-9a-zA-Z+*._> -]+"
type="[()a-zA-Z0-9_ *]+"
if false ;then #use mycalloc, mymalloc, myrealloc to avoid memory cask errors.
	#reformat calloc(s, sizeof(b)) by mycalloc(s, b)
	sed -i "" -E "s|(calloc\($var),[ ]*sizeof\(($type)\)\);|my\1,\2\);|g" $@
	#reformat malloc(s*sizeof(b) by mymalloc(s, b)
	sed -i "" -E "s|(malloc\($var)\*[ ]*sizeof\(($type)\)\);|my\1,\2);|g" $@
	#reformat malloc(sizeof(b)*s by mymalloc(s, b)
	sed -i "" -E "s|(malloc)\(sizeof\(($type)\)[ ]*\*[ ]*($var)\);|my\1(\3,\2);|g" $@
	#reformat realloc(p, s*size(b)) by myrealloc(p, s, b)
	sed -i "" -E "s|(realloc\($var),($var)\*[ ]*sizeof\(($type)\)[ ]*\);|my\1,\2,\3);|g" $@
	#reformat realloc(p, sizeof(b)*s by myrealloc(s, b)
	sed -i "" -E "s|(realloc\($var),[ ]*sizeof\(($type)\)[ ]*\*[ ]*($var)[ ]*\);|my\1,\3,\2);|g" $@
	#reformat realloc(p, s) by myrealloc(s, b, char)
	sed -i "" -E "s|(realloc\($var),($var)\);|my\1,\2,char);|g" $@
	#delete cast
	sed -i -E "s/(\($type\))(mymalloc|myrealloc|mycalloc)/\2/g" $@
fi

if false ;then
	var="[0-9a-zA-Z_.]+"
	arg="[0-9a-zA-Z_+-]+"
	sed -i -E "s|($var)(->p\[)($arg)(\])|P(\1,\3)|g" $@ #repalce  a->p[i] by P(a,i)
	sed -i -E "s/($var|$var\[$var\])P\(\./P(\1./g" $@ #Fix A[a]->P(.b,i) by P(A[a].b,i)
	sed -i -E "s/($var|$var\[$var\])P\(\./P(\1./g" $@ #repeat

	#sed -i -E "s/([^&]+)&P\(/\1PP(/g" $@ #Fix &P(a,i) by PP(A,i) # Manually fix PP()->p[] by &P()->p[]
	#sed -i -E "s/PP(\($var,$var\)->)/&P\1/g" $@ #Fix PP(A,i)-> by &P(A,i)-> 
	#sed -i -E "s/PP(\($var->$var,$var\)->)/&P\1/g" $@ #Fix PP(A->b,i)-> by &P(A->b,i)-> 

	vars="P\($var,$var\)|P\($var->$var,$var\)|P\($var[$var]\.$var,$var\)|P\($var[$var]\.$var->$var,$var\)|$var->$var|$var\[$var\]|$var\[$var\]\.$var|$var->$var\.$var|$var->$var\[$var\]\.$var"
	#sed -i -E "s|(P\($var,$var\))(->p\[)($arg)(\])|P(\1,\3)|g" $@ #replace P(a,i)->p[j] by P(P(a,i),j)
	#sed -i -E "s|(P\($var->$var,$var\))(->p\[)($arg)(\])|P(\1,\3)|g" $@ #replace P(a->b,i)->p[j] by P(P(a->b,i),j)
	#sed -i -E "s|(P\($var->$var->$var,$var\))(->p\[)($arg)(\])|P(\1,\3)|g" $@ #replace P(a->b,i)->p[j] by P(P(a->b,i),j)
	#sed -i -E "s|($var\[$var\])(->p\[)($arg)(\])|P(\1,\3)|g" $@ #replace a[i]->p[j] by P(a[i],j)
	sed -i -E "s/($vars)(->p\[)($arg)(\])/P(\1,\3)/g" $@ #
	sed -i -E "s/($vars)(->p\[)($vars)(\])/P(\1,\3)/g" $@ #
	sed -i -E "s/($vars)(->p\[)($vars)($arg)(\])/P(\1,\3\4)/g" $@ #

	sed -i -E "s|($var)->P\(|P(\1->|g" $@ #fix A->P(a,i) by P(A->a,i)
	sed -i -E "s|($var)->P\(|P(\1->|g" $@ #repeat
fi
#grep 'P([^)]*,[^)]*,[^)]*)' */*.c
#replace P()->p by P(P())
#sed -i -E "s|(P\([^)]*\))->p([^[a-zA-Z])|P(\1)\2|g" $@
#sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;]+)->p([^[a-zA-Z+-])|P(\1)\2|g" $@

#sed -i -E "s|P\(([[:blank:]]+)|\1P(|g" $@

#replace A->nx by NX(A)
#sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<>%-]+)->nx([^a-zA-Z+-])|NX(\1)\2|g" $@
#replace A->NX(B) by NX(A->B)
#sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<%-]+)->NX\(|NX(\1->|g" $@
#replace A->ny by NY(A)
#sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<>%-]+)->ny([^a-zA-Z+-])|NY(\1)\2|g" $@
#replace A->NY(B) by NY(A->B)
#sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<%-]+)->NY\(|NY(\1->|g" $@

#replace A->base by CELL(A)
sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<>%-]+)->base([^a-zA-Z+-])|CELL(\1)\2|g" $@
#replace A->CELL(B) by CELL(A->B)
sed -i -E "s|([^,\"[:blank:]=+():?/*|!&;<%-]+)->CELL\(|CELL(\1->|g" $@