set xlabel "x"
set ylabel "y"
set grid
set terminal png size 800,600 font "Arial,12"
set output "exp.png"

plot "val.dat" using 1:2 with lines title 'Original',\
    "val.dat" using 1:3 with lines title 'Taylor',\
    "val.dat" using 1:4 with lines title 'Pade'
