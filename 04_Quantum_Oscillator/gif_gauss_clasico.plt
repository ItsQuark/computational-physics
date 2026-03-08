set terminal jpeg size 1280, 720
set output 'comparacion_gauss_clasico_centro.jpg'
set multiplot

# Gráfico 1 en la primera columna
set size 0.5, 0.5
set origin 0, 0.5
set xrange [0:250]
set yrange [0:1] 
plot 'valores_medios_g.dat' u 1:2 w l lw 2 t '<x>_{q}', '' u 1:6 w l lw 2 t 'x_{c}'

# Gráfico 2 en la segunda columna
set size 0.5, 0.5
set origin 0.5, 0.5
set xrange [0:250]
set yrange [-35:35] 
plot 'valores_medios_g.dat' u 1:3 w l lw 2 t '<p>_{q}', '' u 1:7 w l lw 2 t 'p_{c}'

# Gráfico 3 ocupando toda la fila siguiente
set size 1, 0.5
set origin 0, 0
set xrange [0:500]
set yrange [145:150] 
plot 'valores_medios_g.dat' u 1:4 w l lw 2 t '<E>_{q}', '' u 1:8 w l lw 2 t 'E_{c}'


