set term gif animate
set output 'schr1_g.gif'

do for [a=0:500:2] {
    set multiplot layout 3,1

    # Densidad de probabilidad
    set xrange [0:1]
    set yrange [0:10] 
    plot 'dens_prob_g.dat' i a u 1:2 w l lw 2 t 'dens prob'

    # Función de onda (parte real e imaginaria)
    set xrange [0:1]
    set yrange [-3.5:3.5] 
    plot 'dens_prob_g.dat' i a u 1:3 w l lw 2 t 'Re(phi)', '' i a u 1:4 w  l lw 2 t 'Im(phi)'

    # Integral de 0 a L
    set xrange [0:500]
    set yrange [0:2] 
    plot 'normalizado_g.dat' i 0:a u 1:2 w l lw 2 t 'Norma'
}

unset multiplot