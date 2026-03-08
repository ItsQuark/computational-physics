set terminal jpeg size 800, 580
set output 'energia_gauss.jpg'

set xrange [0:500]
set yrange[145:150]
plot 'valores_medios_g.dat' u 1:4 w l lw 2 t '<E>'
