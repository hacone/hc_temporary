set terminal svg

set output "separation.svg"

plot "data.dat" u ($10 - 0.2+0.2*rand(0)-0.1):5 w p pt 7 ps 0.4 title "best", \
     "data.dat" u ($10 + 0.2+0.2*rand(0)-0.1):6 w p pt 7 ps 0.4 title "2nd best"

set output "separation-prom.svg"
plot "data.dat" u ($10 + 0.8*rand(0)-0.4):7 w p pt 7 ps 0.4 title "prom"

set xrange [-0.5:*]
set output "round-variants.svg"
plot "data.dat" u ($10 + 0.8*rand(0)-0.4):9 w p pt 7 ps 0.4 title "#vars"
set xrange [*:*]

set terminal png size 1200,1200
set output "separation.png"

plot "data.dat" u ($10 - 0.2+0.2*rand(0)-0.1):5 w p pt 7 ps 0.4 title "best", \
     "data.dat" u ($10 + 0.2+0.2*rand(0)-0.1):6 w p pt 7 ps 0.4 title "2nd best"

set output "separation-prom.png"
plot "data.dat" u ($10 + 0.8*rand(0)-0.4):7 w p pt 7 ps 0.4 title "prom"

set xrange [-0.5:*]
set output "round-variants.png"
plot "data.dat" u ($10 + 0.8*rand(0)-0.4):($9+0.8*rand(0)-0.4) w p pt 7 ps 0.4 title "#vars"
set xrange [*:*]
