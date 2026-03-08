set term gif animate delay 250
set output '0.01.gif'

set palette grey

set xra[-0.5:99.5]
set yra[-0.5:99.5]

do for [a=1:100] {plot '01.dat' i a matrix with image t 'Evolución T=0.01'}
