set term gif animate delay 3
set output 'interaccion1.gif'

#set label 'h=5' at 0.8, 0.8
set title "Interacción con la Luna"

set style line 1 lc rgb 'blue' lw 2

do for [a=1:865:2] {
    
    set xra [-1.1:1.1]
    set yra [-1.1:1.1]
    set size square 

    plot 'interaccion1.dat' i 0:a u 3:4 w l ls 1 t 'Nave', '' i 0:a u 1:2 w l lw 3 t 'Luna'

}