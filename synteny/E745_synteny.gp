set terminal png tiny size 800,800
set output "E745_synteny.png"
set title "Synteny: E745 vs Aus0004"
set ytics ( \
 "*tig00000002" 1.0, \
 "*tig00000009" 199382.0, \
 "tig00000001" 224329.0, \
 "tig00000012" 2986803.0, \
 "tig00000005" 2991027.0, \
 "tig00000007" 3031039.0, \
 "tig00000004" 3047082.0, \
 "tig00000003" 3061815.0, \
 "tig00000008" 3084586.0, \
 "tig00000006" 3100013.0, \
 "" 3115034 \
)
set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "CP003351.2"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set xrange [1:2952485]
set yrange [1:3115034]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "E745_synteny.fplot" title "FWD" w lp ls 1, \
 "E745_synteny.rplot" title "REV" w lp ls 2
