
set terminal svg size 1000, 600
set xrange [-50:2450]
set yrange [0:50]

set output "hg002_ont_def.svg"
plot "< gawk 'NR>2' ./C19.chrX.var.dat" u (200*$1)+($2)+(0.25*$3):5 w imp

exit

set output "chm13_pac_def.svg"
plot "< gawk 'NR>2' ./vf_chm13_pac_def.dat" u (200*$1)+($2)+(0.25*$3):5 w imp

set yrange [-50:50]
set output "chm13_both_def.svg"
plot "< gawk 'NR>2' ./vf_chm13_ont_def.dat" u (200*$1)+($2)+(0.25*$3):5 w imp, \
     "< gawk 'NR>2' ./vf_chm13_pac_def.dat" u (200*$1)+($2)+(0.25*$3):(-1*$5) w imp

exit

set terminal svg size 1200, 300
set yrange [-10:10]
set output "chm13_both_def_zm.svg"
plot "< gawk 'NR>2' ./vf_chm13_ont_def.dat" u (200*$1)+($2)+(0.25*$3):5 w imp, \
     "< gawk 'NR>2' ./vf_chm13_pac_def.dat" u (200*$1)+($2)+(0.25*$3):(-1*$5) w imp
