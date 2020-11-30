#!/bin/bash

# Init global parameters
honpas_gitrev=`git rev-parse --short HEAD 2>/dev/null`
test -z "${honpas_gitrev}" && honpas_gitrev="unknown"

# Init system
sysname=`/bin/ls -t | grep '\.out$' | head -n 1 | sed -e 's/\.out$//'`
bands_src=`/bin/ls -t | grep '\.bands$' | head -n 1`
bands_dst="${sysname}.bands.dat"
bands_gnu="${sysname}.bands.gnu"
bands_pdf="${sysname}.bands.pdf"
bands_ref_src="../Refs/${sysname}.bands"
bands_ref_dst="${sysname}.bands.ref"
bands_title=`echo "System: ${sysname}" | sed -e 's,_,/,g'`
bands_title="${bands_title} - Git revision: ${honpas_gitrev}"

# Convert data
echo "Generating Gnuplot input data"
echo "===> ${bands_ref_dst}"
gnubands -o "${bands_ref_dst}" "${bands_ref_src}"
echo "===> ${bands_dst}"
gnubands -o "${bands_dst}" "${bands_src}"

# Find data range
xmax=`awk '{print $1}' "${bands_dst}" | sed -e '/^#.*/d' | sort | tail -n 1`

# Generate Gnuplot script
cat >"${bands_gnu}" <<EOF
set terminal pdf color
set output "${bands_pdf}"

set title "{${bands_title}}"
set xrange [0.0:${xmax}]
#set yrange [-20.0:15.0]
set xtics ("L"  0.000000, "Γ"  0.530288, "X"  1.142611, "W"  1.448773, "K"  1.665262, "Γ"  2.314730)

plot "${bands_dst}" using 1:2 title "Now" with lines, \
     "${bands_ref_dst}" using 1:2 title "秦新明" with lines dashtype 2
EOF

# Generate PDF plot
gnuplot "${bands_gnu}"
