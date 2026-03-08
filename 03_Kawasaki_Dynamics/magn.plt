# m(T) para N = 32

set terminal jpeg
set output 'magn_32.jpg'

m_min = -0.6
m_max = 0.6
T = 7

set xrange [0:T]
set yrange [m_min:m_max]
set xlabel "T"
set ylabel "m(T)"
set xtics 1
set ytics 0.1
set key
set grid
set title "N=32"

plot '32.dat' u 1:2 w l t '<m>_{sup, N=32}', '' u 1:3 w l t '<m>_{inf, N=32}', \
     # '64.dat' u 1:2 w l t '<m>_{sup, N=64}', '' u 1:3 w l t '<m>_{inf, N=64}', \
     # '128.dat' u 1:2 w l t '<m>_{sup, N=128}', '' u 1:3 w l t '<m>_{inf, N=128}'

unset output


# m(T) para N = 64

set terminal jpeg
set output 'magn_64.jpg'
set title "N = 64"

plot '64.dat' u 1:2 w l t '<m>_{sup, N=64}', '' u 1:3 w l t '<m>_{inf, N=64}'
unset output


# m(T) para N = 128

set terminal jpeg
set output 'magn_128.jpg'
set title "N = 128"

plot '128.dat' u 1:2 w l t '<m>_{sup, N=128}', '' u 1:3 w l t '<m>_{inf, N=128}'
unset output



# IDEM PARA M=N^2/2

# m(T) para N = 32

set terminal jpeg
set output 'magn_32m.jpg'

m_min = -0.3
m_max = 0.8
T = 20

set xrange [0:T]
set yrange [m_min:m_max]
set xlabel "T"
set ylabel "m(T)"
set xtics 2
set ytics 0.1
set key
set grid
set title "N=32"

plot '32m.dat' u 1:2 w l t '<m>_{sup, N=32}', '' u 1:3 w l t '<m>_{inf, N=32}', \
     # '64.dat' u 1:2 w l t '<m>_{sup, N=64}', '' u 1:3 w l t '<m>_{inf, N=64}', \
     # '128.dat' u 1:2 w l t '<m>_{sup, N=128}', '' u 1:3 w l t '<m>_{inf, N=128}'

unset output


# m(T) para N = 64

set terminal jpeg
set output 'magn_64m.jpg'
set title "N = 64"

plot '64m.dat' u 1:2 w l t '<m>_{sup, N=64}', '' u 1:3 w l t '<m>_{inf, N=64}'
unset output


# m(T) para N = 128

set terminal jpeg
set output 'magn_128m.jpg'
set title "N = 128"

plot '128m.dat' u 1:2 w l t '<m>_{sup, N=128}', '' u 1:3 w l t '<m>_{inf, N=128}'
unset output
