set term gif animate
set output 'prob.gif'

do for [a=0:100:2] {
    set multiplot layout 3,1

    # Densidad de probabilidad n=0
    set xrange [0:1]
    set yrange [0:7] 
    set arrow 1 from 0.4, -0.3 to 0.4, 0.3 nohead lw 3
    set arrow 2 from 0.6, -0.3 to 0.6, 0.3 nohead lw 3
    set label 'n=0 ; ω=200' at 0.03, 5.5
    plot 'prob_n0.dat' i a u 1:2 w l lw 2 t 'P_{n}', '' i a u 1:3 w l lw 2 t 'P_{c}'
    unset label
    unset arrow 1
    unset arrow 2

    # Densidad de probabilidad n=3
    set xrange [0:1]
    set yrange [0:7]
    set arrow 3 from 0.235, -0.3 to 0.235, 0.3 nohead lw 3
    set arrow 4 from 0.765, -0.3 to 0.765, 0.3 nohead lw 3
    set label 'n=3 ; ω=200' at 0.03, 5.5
    plot 'prob_n3.dat' i a u 1:2 w l lw 2 t 'P_{n}', '' i a u 1:3 w l lw 2 t 'P_{c}'
    unset label
    unset arrow 3
    unset arrow 4

    # Densidad de probabilidad n25
    set xrange [0:1]
    set yrange [0:7] 
    set arrow 5 from 0.118, -0.3 to 0.118, 0.3 nohead lw 3
    set arrow 6 from 0.882, -0.3 to 0.882, 0.3 nohead lw 3
    set label 'n=25 ; ω=700' at 0.03, 5.5
    plot 'prob_n25.dat' i a u 1:2 w l lw 2 t 'P_{n}', '' i a u 1:3 w l lw 2 t 'P_{c}'
    unset label
    unset arrow 5
    unset arrow 6
}
unset multiplot
