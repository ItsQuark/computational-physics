set output 'grafico.png'

plot 'energia2.dat' u 1:2 w lp ps 3 t '<e>(T) exp', '' u 1:3 w lp ps 3 t '<e>(T) teo'