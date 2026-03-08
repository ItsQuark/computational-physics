# Establezco los límites de los ejes x e y, en función de la longitud de red N
N = 128
xmin = 0 - 0.5
xmax = N - 0.5
ymin = 0 - 0.5
ymax = N - 0.5

# Cantidad de fotogramas (pasos montecarlo máximos)
pm = 500

# Primer gif
set term gif animate delay 5
set output '1.gif'

set palette grey
set title "T = 0.01"
set size square
set xra[xmin:xmax]
set yra[ymin:ymax]

do for [a=1:pm] {plot '1.dat' i a matrix with image }
unset output

# Segundo gif
set term gif animate delay 9
set output '2.gif'

set palette grey
set title "T = 3.00"
set size square
set xra[xmin:xmax]
set yra[ymin:ymax]

do for [a=1:pm] {plot '2.dat' i a matrix with image }
unset output

# Tercer gif
set term gif animate delay 5
set output '3.gif'

set palette grey
set title "T = 8.00"
set size square
set xra[xmin:xmax]
set yra[ymin:ymax]

do for [a=1:pm] {plot '3.dat' i a matrix with image }
unset output

# Cuarto gif
set term gif animate delay 5
set output '4.gif'

set palette grey
set title "T = 2.40"
set size square
set xra[xmin:xmax]
set yra[ymin:ymax]

do for [a=1:pm] {plot '4.dat' i a matrix with image }
unset output