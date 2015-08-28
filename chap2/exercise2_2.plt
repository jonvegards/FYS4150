set terminal svg enhanced size 1000 1000 fname "Times" fsize 36
set output "exercise2_2.svg"
set title "Log-log plot of the error"
set xlabel "x"
set ylabel "y"
plot "./exercise2_2.dat" using 1:2 title ""