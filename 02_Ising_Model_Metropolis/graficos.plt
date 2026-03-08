set terminal jpeg
set output 'salida.jpg'

plot 'energia.dat' u 1:2 w lp ps 3

# Luego pon en la terminal lo de siempre: gnuplot graficos.plt
