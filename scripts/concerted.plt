
# for global variation distribution
# plot "log" u (4*200*$1 + 4*$2 + $3):5 pt 3

# set terminal qt font "Sans,14"

set terminal svg font "Sans,14"

set grid
set ylabel "%Div"
set xlabel "Units Separated"

set output "DATA.all.svg"
plot "DATA" u 1:2 w l lw 3 title "On all sites"

set output "DATA.snvs.svg"
plot "DATA" u 1:3 w l lw 3 title "On SNV sites"
