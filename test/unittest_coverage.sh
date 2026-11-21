#!/bin/sh
#Print the line coverage of sys/, math/, and lib/ by the unit tests as a markdown
#table. Requires a build configured with --enable-gcov, and make check to have
#been run. Only gcov is used. Called by unittest_report.sh.
#Usage: unittest_coverage.sh TOP_SRCDIR TOP_BUILDDIR
top_srcdir=$(cd "$1" && pwd -P)
top_builddir=$(cd "$2" && pwd -P)
dirs="sys math lib"

if ! command -v gcov >/dev/null 2>&1 ;then
	echo "Not available: gcov is not found."
	exit 0
fi
#Collect "source|lines executed|lines" for every object file with coverage data.
#A template source (e.g. math/mat.c) is included by several objects; the
#object with the most lines executed is used.
list=$(for dir in $dirs ;do
	for gcda in "$top_builddir/$dir"/*.gcda "$top_builddir/$dir"/.libs/*.gcda ;do
		[ -f "$gcda" ] || continue
		objdir=$(dirname "$gcda")
		(cd "$top_builddir/$dir" && gcov -n -o "$objdir" "$gcda" 2>/dev/null) | awk -v dir="$top_builddir/$dir" '
			/^File /{fn=$0; sub(/^File \047/, "", fn); sub(/\047$/, "", fn); if(fn !~ /^\//) fn=dir "/" fn; next}
			/^Lines executed:/{
				if(fn=="") next
				split($0, a, /[:% ]+/) #Lines executed:85.00% of 120
				pct=a[3]; n=a[5]
				printf "%s|%d|%d\n", fn, int(pct*n/100+0.5), n
				fn=""
			}'
	done
done)
if [ -z "$list" ];then
	echo "Not available: configure with \`--enable-gcov\` and run \`make check\` or \`make doc\` to measure coverage."
	exit 0
fi
#Normalize paths, keep files in the source tree, and take the maximum per file.
printf '%s\n' "$list" | while IFS='|' read -r fn hit n ;do
	d=$(cd "$(dirname "$fn")" 2>/dev/null && pwd -P) || continue
	case "$d/" in
		"$top_srcdir"/*) echo "${d#$top_srcdir/}/$(basename "$fn")|$hit|$n" ;;
	esac
done | awk -F'|' -v dirs="$dirs" '
	$1 ~ /\.c$/ {
		split($1, p, "/"); if(index(" " dirs " ", " " p[1] " ")==0) next
		if(!($1 in tot) || $2>hit[$1]){hit[$1]=$2; tot[$1]=$3}
	}
	END{
		ndir=split(dirs, d, " ")
		for(f in tot){split(f, p, "/"); dh[p[1]]+=hit[f]; dt[p[1]]+=tot[f]; th+=hit[f]; tt+=tot[f]}
		printf "Line coverage of %s by the unit tests (gcov): %.1f%% (%d of %d lines).\n\n", dirs, tt?100*th/tt:0, th, tt
		printf "| Folder | Coverage | Lines executed | Lines |\n|--------|----------|----------------|-------|\n"
		for(i=1; i<=ndir; i++) if(dt[d[i]]) printf "| `%s` | %.1f%% | %d | %d |\n", d[i], 100*dh[d[i]]/dt[d[i]], dh[d[i]], dt[d[i]]
		printf "\n| File | Coverage | Lines executed | Lines |\n|------|----------|----------------|-------|\n"
		n=0; for(f in tot) names[++n]=f
		for(i=2; i<=n; i++){x=names[i]; for(j=i-1; j>0 && names[j]>x; j--) names[j+1]=names[j]; names[j+1]=x}
		for(i=1; i<=n; i++){f=names[i]; if(tot[f]) printf "| `%s` | %.1f%% | %d | %d |\n", f, 100*hit[f]/tot[f], hit[f], tot[f]}
	}'
