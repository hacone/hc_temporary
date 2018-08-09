
set terminal png

set xrange [0:22.5]
set yrange [0:30]
set output "vf-gap.png"
plot \
"< gawk '$7 == 0' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 0 title "vf - gap 0",\
"< gawk '$7 == 1' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 1 title "1"

#,\
#"< gawk '$7 == 2' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 2 title "2",\
#"< gawk '$7 == 3' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 3 title "3",\
#"< gawk '$7 == 4' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 4 title "4",\
#"< gawk '$7 == 5' ./.stats" u 5:4 with points ps 0.7 pt 7 lc 5 title "5" 

set xrange [0:22.5]
set yrange [45:90]
set output "vf-scr.png"
plot \
"< gawk '$7 == 0' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 0 title "vf - scr 0",\
"< gawk '$7 == 1' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 1 title "1"
#,\
#"< gawk '$7 == 2' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 2 title "2",\
#"< gawk '$7 == 3' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 3 title "3",\
#"< gawk '$7 == 4' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 4 title "4",\
#"< gawk '$7 == 5' ./.stats" u 5:3 with points ps 0.7 pt 7 lc 5 title "5" 

set xrange [0:15]
set yrange [45:90]
set output "ovlp-scr.png"
plot \
"< gawk '$7 == 0' ./.stats" u ($6+(rand(0)-0.5)*0.4):3 with points ps 0.7 pt 7 lc 0 title "ovlp - scr 0",\
"< gawk '$7 == 1' ./.stats" u ($6+(rand(0)-0.5)*0.4):3 with points ps 0.7 pt 7 lc 1 title "1"

#,\
#"< gawk '$7 == 2' ./.stats" u 6:3 with points ps 0.7 pt 7 lc 2 title "2",\
#"< gawk '$7 == 3' ./.stats" u 6:3 with points ps 0.7 pt 7 lc 3 title "3",\
#"< gawk '$7 == 4' ./.stats" u 6:3 with points ps 0.7 pt 7 lc 4 title "4",\
#"< gawk '$7 == 5' ./.stats" u 6:3 with points ps 0.7 pt 7 lc 5 title "5" 


set xrange [45:90]
set yrange [0:30]
set output "scr-gap.png"
plot \
"< gawk '$7 == 0' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 0 title "scr - gap 0",\
"< gawk '$7 == 1' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 1 title "1"
#,\
#"< gawk '$7 == 2' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 2 title "2",\
#"< gawk '$7 == 3' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 3 title "3",\
#"< gawk '$7 == 4' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 4 title "4",\
#"< gawk '$7 == 5' ./.stats" u 3:4 with points ps 0.7 pt 7 lc 5 title "5" 

