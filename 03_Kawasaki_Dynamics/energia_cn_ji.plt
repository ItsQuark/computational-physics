# <e>(T) para N={32, 64, 128}

T_max = 7

set terminal jpeg 
set output 'energia.jpg'
set xrange [0:T_max]
set xlabel "T"
set ylabel "<e>(T)"
set xtics 2
set grid
set key at screen 0.94, 0.22



plot '32.dat' u 1:4 w l t '<e>_{N=32}', \
     '64.dat' u 1:4 w l t '<e>_{N=64}', \
     '128.dat' u 1:4 w l t '<e>_{N=128}'

unset output


# cn(T) para N={32, 64, 128}

set terminal jpeg 
set output 'cn.jpg'
set ylabel "c_n(T)"
set yrange [0:4.3]
set key at screen 0.94, 0.95

plot '32.dat' u 1:5 w l t 'c_{N=32}', \
     '64.dat' u 1:5 w l t 'c_{N=64}', \
     '128.dat' u 1:5 w l t 'c_{N=128}'

unset output


# ji(T) para N={32, 64, 128}

set terminal jpeg 
set output 'ji.jpg'
set ylabel "χ_n(T)"
set yrange [0:3]
set key at screen 0.94, 0.95

plot '32.dat' u 1:6 w l t 'χ_{N=32}', \
     '64.dat' u 1:6 w l t 'χ_{N=64}', \
     '128.dat' u 1:6 w l t 'χ_{N=128}'

unset output

# IDEM PARA M=N^2/2

# <e>(T) para N={32, 64, 128}

T_max = 20

set terminal jpeg 
set output 'energia_m.jpg'
set xrange [0:T_max]
set yrange [-2:-0.5]
set xlabel "T"
set ylabel "<e>(T)"
set xtics 2
set grid
set key at screen 0.94, 0.22



plot '32m.dat' u 1:4 w l t '<e>_{N=32}', \
     '64m.dat' u 1:4 w l t '<e>_{N=64}', \
     '128m.dat' u 1:4 w l t '<e>_{N=128}'

unset output


# cn(T) para N={32, 64, 128}

set terminal jpeg 
set output 'cn_m.jpg'
set ylabel "c_n(T)"
set yrange [0:2.5]
set key at screen 0.94, 0.95

plot '32m.dat' u 1:5 w l t 'c_{N=32}', \
     '64m.dat' u 1:5 w l t 'c_{N=64}', \
     '128m.dat' u 1:5 w l t 'c_{N=128}'

unset output


# ji(T) para N={32, 64, 128}

set terminal jpeg 
set output 'ji_m.jpg'
set ylabel "χ_n(T)"
set yrange [0:1.6]
set key at screen 0.94, 0.95

plot '32m.dat' u 1:6 w l t 'χ_{N=32}', \
     '64m.dat' u 1:6 w l t 'χ_{N=64}', \
     '128m.dat' u 1:6 w l t 'χ_{N=128}'

unset output