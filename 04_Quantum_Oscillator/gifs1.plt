set term gif animate delay 5
set output 'med_n3.gif'

do for [a=0:500:2] {
    set multiplot layout 3,1

    # Valor medio de x
    set xrange [0:500]
    set yrange [0:1] 
    set label 1 'n=3' at 20, 0.8
    plot 'valores_medios.dat' i 0:a u 1:2 w l lw 2 t '<x>'

    # Valores medios de p y H
    set xrange [0:500]
    set yrange [-5:1400] 
    unset label
    plot 'valores_medios.dat' i 0:a u 1:3 w l lw 2 t '<p>', '' i 0:a u 1:4 w  l lw 2 t '<H>', '' i 0:a u 1:5 w  l lw 2 t '<H>_{teo}' 

    # Δx·Δp
    set xrange [0:500]
    set yrange [0:7] 
    unset label
    plot 'valores_medios.dat' i 0:a u 1:6 w l lw 2 t 'Δx·Δp', '' i 0:a u 1:7 w l lw 2 t '(Δx·Δp)_{teo}'
    unset label
}

unset multiplot