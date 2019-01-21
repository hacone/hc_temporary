
set terminal svg size 1200,400

set output "xxxxx.svg"

plot "xxxxx" u (2*(4*200*$1 + 4*$2 + $3)):5 w imp

#pt 3 ps 0.8
