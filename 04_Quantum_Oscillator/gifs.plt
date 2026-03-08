set term gif animate delay 5
set output 'dens_n2.gif'

do for [a=0:500:2] {
    set multiplot layout 3,1

    # Densidad de probabilidad
    set xrange [0:1]
    set yrange [0:8] 
    set label 'n=2' at 0.03, 6.5
    plot 'dens_prob.dat' i a u 1:2 w l lw 2 t 'dens prob'

    # Función de onda (parte real e imaginaria)
    set xrange [0:1]
    set yrange [-3.5:3.5] 
    unset label
    plot 'dens_prob.dat' i a u 1:3 w l lw 2 t 'Re(phi)', '' i a u 1:4 w  l lw 2 t 'Im(phi)'

    # Integral de 0 a L
    set xrange [0:500]
    set yrange [0:2] 
    unset label
    plot 'normalizado.dat' i 0:a u 1:2 w l lw 2 t 'Norma'
}

unset multiplot