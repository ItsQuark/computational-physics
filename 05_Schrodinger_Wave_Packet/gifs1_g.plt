set term gif animate delay 5
set output 'valores_medios_g.gif'

do for [a=0:500:2] {
    set multiplot layout 3,1

    # Valor medio de x
    set xrange [0:500]
    set yrange [0:1] 
    plot 'valores_medios_g.dat' i 0:a u 1:2 w l lw 2 t '<x>'

    # Valores medios de p y H
    set xrange [0:500]
    set yrange [-50:50] 
    plot 'valores_medios_g.dat' i 0:a u 1:3 w l lw 2 t '<p>', '' i 0:a u 1:4 w  l lw 2 t '<H>'

    # Integral de 0 a L
    set xrange [0:500]
    set yrange [0:20]
    unset label
    plot 'valores_medios_g.dat' i 0:a u 1:5 w l lw 2 t 'Δx·Δp'
}

unset multiplot