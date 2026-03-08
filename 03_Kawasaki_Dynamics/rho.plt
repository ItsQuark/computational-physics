# ρ(y) para N = 32

set terminal jpeg
set output '32_rhoM.jpg'

set xrange [0:32]
set yrange [-1.1:1.1]
set xlabel "y"
set ylabel "ρ(y)"
set ytics 0.25
set key at screen 0.94, 0.25
set grid
set title "N = 32"

plot '32_1rho.dat' w l t 'T = 0.01', \
     '32_2rho.dat' w l t 'T = 1.50', \
     '32_3rho.dat' w l t 'T = 2.00', \
     '32_4rho.dat' w l t 'T = 8.00'
unset output


# ρ(y) para N = 64

set terminal jpeg
set output '64_rhoM.jpg'

set xrange [0:64]
set yrange [-1.1:1.1]
set xlabel "y"
set ylabel "ρ(y)"
set key at screen 0.94, 0.25
set ytics 0.25
set grid
set title "N = 64"

plot '64_1rho.dat' w l t 'T = 0.01', \
     '64_2rho.dat' w l t 'T = 1.50', \
     '64_3rho.dat' w l t 'T = 2.00', \
     '64_4rho.dat' w l t 'T = 8.00'
unset output


# ρ(y) para N = 128

set terminal jpeg
set output '128_rhoM.jpg'

set xrange [0:128]
set yrange [-1.1:1.1]
set xlabel "y"
set ylabel "ρ(y)"
set key at screen 0.94, 0.25
set ytics 0.25
set grid
set title "N = 128"

plot '128_1rho.dat' w l t 'T = 0.01', \
     '128_2rho.dat' w l t 'T = 1.50', \
     '128_3rho.dat' w l t 'T = 2.00', \
     '128_4rho.dat' w l t 'T = 8.00'
unset output
