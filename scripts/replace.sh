var="[()0-9a-zA-Z+*._> -]+"
type="[()a-zA-Z0-9_ *]+"
#minus sign must be the last or first inside a []. no need to escape ><
#reformat calloc(s, sizeof(b)) by mycalloc(s, b)
sed -i "" -E "s|(calloc\($var),[ ]*sizeof\(($type)\)\);|my\1,\2\);|g" $@
#reformat malloc(s*sizeof(b) by mymalloc(s, b)
sed -i "" -E "s|(malloc\($var)\*[ ]*sizeof\(($type)\)\);|my\1,\2);|g" $@
#reformat malloc(sizeof(b)*s by mymalloc(s, b)
sed -i "" -E "s|(malloc)\(sizeof\(($type)\)[ ]*\*[ ]*($var)\);|my\1(\3,\2);|g" $@
#reformat realloc(p, s*size(b)) by myrealloc(p, s, b)
sed -i "" -E "s|(realloc\($var),($var)\*[ ]*sizeof\(($type)\)[ ]*\);|my\1,\2,\3);|g" $@
#reformat remalloc(p, sizeof(b)*s by mymalloc(s, b)
sed -i "" -E "s|(realloc\($var),[ ]*sizeof\(($type)\)[ ]*\*[ ]*($var)[ ]*\);|my\1,\3,\2);|g" $@
#reformat remalloc(p, s by mymalloc(s, b, char)
sed -i "" -E "s|(realloc\($var),($var)\);|my\1,\2,char);|g" $@
