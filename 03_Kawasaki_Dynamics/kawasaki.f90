program kawasaki
    use randomnumber
    implicit none
    integer :: n, dim_T, w, v, magn, i
    integer*8 :: fin_pasos, archivo_append
    real*8 :: T_min, T_max, max_cn, max_ji, Tc_cn, Tc_ji
    integer, allocatable :: s(:,:)
    real*8, allocatable :: T(:), Tc(:)
    character(len=1) :: eleccion


    ! Distingo la funcionalidad del ejecutable: generación de gifs (fuerzo pocos pasos MC) o cálculo de promedios (doy la posibilidad de dar más pasos MC)
    
    print*, "¿Desea escribir los datos en un fichero necesario para crear los gifs (f) o calcular los promedios (p)?"
    ! f: por defecto, el algoritmo usa 100 pasos montecarlo debido al peso de los archivos (100 pasos MC = 18.76 MB)
    ! Actualización: ahora escribo cada 10 pasos MC, por lo que poner fin_pasos = 1000 sigue equivaliendo a 18.76 MB 
    fin_pasos = 5000
    ! p: si elige promediar, no se escribe nada en un fichero para no saturarlos, y doy la posibilidad de elegir muchos pasos montecarlo
    read(*,*) eleccion
    if (eleccion == "p") then
        print*, "Pasos montecarlo (numero entero par mayor que cero):"
        read(*,*) fin_pasos
    end if

    call dran_ini(1968339)

    if ( eleccion=="f" ) then

        ! Si solo se desea hacer los gifs, ejecutamos esta parte del código en la que la longitud de red n permanece fija
        ! Longitud de la red (# puntos = n²) (solo funciona para n par, de modo que se permita una magnetización nula)
        n = 128
        i = 3
        !magn = 0
        magn = n*n/2

        ! Creación del vector que almacena las temperaturas. Para la elección 'f', son 4 cuatro temperaturas, y 3 las T críticas, designadas por i (N=32 -> i=1, N=64 -> i=2, N=128 -> i=3)
        dim_T = 4
        allocate (T(dim_T), Tc(3))
        
        !Tc = [2.6349494949494949d0, 2.6349494949494949d0, 2.8368686868686863d0]
        !T = [0.6*Tc(i), 0.7*Tc(i), 0.8*Tc(i), Tc(i)]                           ! Tc incierta
        T = [0.01d0, 3d0, 8d0, 2.4d0]

        allocate (s(n,n))
        do w = 1, dim_T

            ! Llamada a la subrutina que inicializa la matriz con la magnetización M deseada. Debe ser par, de otro modo no existe el número de espines positivos tal que de esa magnetización impar
            ! La magnetización puede tomar los valores desde -n*n+2*n a n*n-2*n, siendo par, para asegurar que haya el suficiente número de espines para llenar las filas 1 y n
            call init(magn)
    
            ! do v = 1, n
            !     print*, s(v, :)
            ! end do
            ! print*, ""
            ! print*, ""

            ! Llamada al algoritmo metropolis (input: índice w que etiqueta la temperatura deseada) tantas veces como temperaturas haya
            call metropolis(w)
        end do
    
    else if (eleccion=="p") then

        ! Si se desea hacer cálculo se ejecuta esta parte, donde varían tanto la temperatura como n=[32, 64, 128] a M/(n*n) = m fijo

        n=32             ! No tocar
        !magn = n*n/2     ! Si M distinto de 0, decomentar
        magn = 0        ! Si M = 0, decomentar

        ! Creador del vector T: cuánto mayor sea dim_T, más datos de M(T) (y por tanto más tarda el finalizar el programa)
        ! Para hacer apartado 2, usa2r dim_T grande (>10)
        ! Para hacer apartado 3, usar dim_T=4
        dim_T = 20
        T_min = 2.00d0
        T_max = 3.00d0
        allocate(T(dim_T))

        do w = 1, dim_T
            T(w) = T_min + (w-1)*(T_max-T_min)/(dim_T-1)
        end do
        
        ! M=0:

        do v = 1, 3

            ! Para encontrar el máximo en cn(T) y ji(T)
            max_cn = 0d0
            max_ji = 0d0

            allocate (s(n,n)) 
            do w = 1, dim_T
                ! IMPORTANTE, EL ARCHIVO QUE SE ABRE PARA CADA N DEBE SER ABIERTO EN MODO APPEND CUANDO SE ABRE PARA LA SEGUNDA TEMPERATURA EN ADELANTE (w>1)
                archivo_append = 1
                if (w==1) archivo_append = 0
                call init(magn)
                call metropolis(w)
                print*, "Temperatura", T(w), "hecho"
            end do
            print*, "Longitud de red", n, "hecho"
            
            write(*,*) n, Tc_cn, Tc_ji

            n = n+n
            !magn = n*n/2
            magn = 0
            deallocate(s)

        end do

        ! M != 0:

        n=32             ! No tocar
        magn = n*n/2

        do v = 1, 3

            ! Para encontrar el máximo en cn(T) y ji(T)
            max_cn = 0d0
            max_ji = 0d0

            allocate (s(n,n)) 
            do w = 1, dim_T
                ! IMPORTANTE, EL ARCHIVO QUE SE ABRE PARA CADA N DEBE SER ABIERTO EN MODO APPEND CUANDO SE ABRE PARA LA SEGUNDA TEMPERATURA EN ADELANTE (w>1)
                archivo_append = 1
                if (w==1) archivo_append = 0
                call init(magn)
                call metropolis(w)
            end do
            
            write(*,*) n, Tc_cn

            n = n+n
            magn = n*n/2
            deallocate(s)

        end do


    else 
        print*, "Solo puede insertar 'f' o 'p'"
    end if
    
    stop
    contains

    ! Subrutina init: inicializa la red del modelo de Ising con la magnetización M deseada, con spin -1 en la fila 1, y +1 en la fila N

    subroutine init(magnetizacion)
        implicit none
        integer :: k, l, magnetizacion, espines_positivos
        ! integer:: desorden, pos11, pos12, pos21, pos22, aux

        ! Sabiendo M, conocemos el número de espines positivos mediante: S(+) = 1/2 * (M + N*N)
        espines_positivos = (magnetizacion+n*n)/2

        if ( espines_positivos<2*n-n*n ) then
            print*, "Error, debe introducir una magnetización mayor o igual a 2*n-n*n"
            stop
        end if

        ! Primero lleno la matriz de espines negativos
        s=-1

        ! Coloco +1 'espines_positivos' veces de forma ordenada
        do k = n, 1, -1
            do l = n, 1, -1
                if ( espines_positivos>0 ) then
                    s(k,l) = 1
                    espines_positivos = espines_positivos-1
                else
                    exit
                end if
            end do
        end do

        !Mezclo las orientaciones de los espines intercambiando posiciones de forma aleatoria, 'desorden' veces
        ! desorden = 10000
        ! do k = 1, desorden
        !     pos11 = i_dran(int(n-2, 8))+1     ! Para que no desordene las filas 1 y N
        !     pos12 = i_dran(int(n, 8))
        !     pos21 = i_dran(int(n-2, 8))+1
        !     pos22 = i_dran(int(n, 8))

        !     aux = s(pos11, pos12)
        !     s(pos11, pos12) = s(pos21, pos22)
        !     s(pos21, pos22) = aux
        ! end do

    end subroutine

    ! Subrutina metropolis. Input: indice_T ( se usa para etiquetar las temperaturas en T(indice_T) )

    subroutine metropolis(indice_T)
        implicit none
        real*8 :: p
        integer*8 :: delta_E, iter, pasos_MC, contador
        integer :: indice_T, unit_number, k, l, vecino, k_prima, l_prima, aux, unit_number2
        character(len=20) :: filename, filename2

        integer :: max_fils_neg
        integer*8 :: magn_inf, magn_sup
        real*8 :: promedio_magn_inf, promedio_magn_sup, promedio_magn_sup_2
        real*8 :: E, promedio_E, promedio_EE, cn, ji

        ! Ficheros donde almaceno la evolucion de la red de spines para distintas temperaturas, y ciertos datos para distintas n. Los ficheros se crean dinámicamente.

        unit_number = 10+indice_T                   ! indice_T llega a dim_T=4 como mucho (unit_number = 11, 12, 13, 14)
        unit_number2 = 10+dim_T+v                   ! dim_T=4 y v etiqueta a n (unit_number2 = 15, 16, 17), de modo que: v=1 -> 32      v=2 -> 64       v=3 -> 128   
        write(filename, *) indice_T
        write(filename2, *) n
        filename = trim(adjustl(filename)) // '.dat'
        
        if ( magn==0 ) then
            filename2 = trim(adjustl(filename2)) // '.dat'
        else
            filename2 = trim(adjustl(filename2)) // 'm.dat'
        end if

        ! Ficheros donde se almacenan para distintas temperatuas la evolución de la red de espines (1.dat, 2.dat, ...)
        if (eleccion=='f') open(unit=unit_number, file=filename, status="unknown")
        
        ! Ficheros donde se almacenan para cada n fija, variando T: m_
        if (eleccion=='p') then
            if (archivo_append == 0) then
                open(unit=unit_number2, file=filename2, status="unknown")
            else if (archivo_append == 1) then
                open(unit=unit_number2, file=filename2, status="unknown", position="append")   
            end if        
        end if
        
        ! Antes de iniciar el bucle que itera sobre pasos MC, se inicializan las variables que almacenan promedios

        contador = 0
        promedio_magn_inf = 0d0
        promedio_magn_sup = 0d0
        promedio_E = 0d0
        promedio_EE = 0d0
        promedio_magn_sup_2 = 0d0
        
        do pasos_MC = 0, fin_pasos-1

            ! Escribo la evolucion de la matriz para la temperatura deseada en el fichero correspondiente, si el usuario eligió 'f'

            if (eleccion == 'f') then
                if ( mod(pasos_MC, 10)==0 ) then
                    do k=n, 1, -1
                        ! La escribo al revés porque por razones que desconozco Gnuplot representa los datos 'del revés'
                        write(unit_number,*) s(k,:)
                    end do
                    write(unit_number,*)
                    write(unit_number,*)
                end if
            end if


            ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ! EVOLUCIÓN: Se hacen n*n preguntas en casillas aleatorias para que la red evolucione 1 paso MC 
            ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

            do iter = 1, n*n

                ! Elijo un punto al azar de la red (que no pertenezca a las filas 1 y n para no alterar los espines fijos)
                k = i_dran(int(n-2,8))+1
                l = i_dran(int(n,8))

                ! Elijo el punto vecino al anterior:    1 -> arriba         2 -> derecha        3 -> izquierda          4 -> abajo
                ! Se impone que el spin vecino no puede caer en las filas 1 y n, y las condiciones de contorno periódicas sobre el eje x

                if ( (k.ne.2).AND.(k.ne.(n-1)) ) then                     
                    vecino = i_dran(int(4, 8))
                else if ( k==2 ) then
                    vecino = i_dran(int(3, 8))+1
                else
                    vecino = i_dran(int(3, 8))
                end if

                if ( vecino==1 ) then
                    k_prima = k-1
                    l_prima = l
                else if ( vecino==2 ) then
                    k_prima = k
                    if ( l==n ) then
                        l_prima = 1
                    else
                        l_prima = l+1
                    end if
                else if ( vecino==3 ) then
                    k_prima = k
                    if ( l==1 ) then
                        l_prima = n
                    else
                        l_prima = l-1
                    end if
                else if ( vecino==4 ) then
                    k_prima = k+1
                    l_prima = l
                end if

                ! Una vez tengo dos puntos vecinos seleccionados mediante (k,l) y (k_prima,l_prima), compruebo si tienen el mismo valor. Si es así, delta_E=0 y no
                ! es necesario discutir si hay intercambio de valores o no. Si no, calculo delta_E y evaluo p.

                if ( s(k,l)==s(k_prima, l_prima) ) then
                    delta_E = 0d0
                else
                    delta_E = g(k, l, k_prima, l_prima)

                    ! Se determina la probabilidad de transición p
                    p = min(1.d0, exp(-delta_E/T(indice_T)))
                    
                    ! Si un número aleatorio con pdf uniforme generado entre 0 y 1 sale estrictamente menor que p, se intercambian los espines
                    if ( dran_u()<p ) then
                        aux = s(k,l)
                        s(k,l) = s(k_prima, l_prima)
                        s(k_prima, l_prima) = aux
                    end if

                end if

                

            end do 
            ! Finaliza un paso MC
            
            ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
            ! CÁLCULOS: se calcula la energía y la magnetización, cada 5 pasos MC y a partir del momento en el que la evolución va por la mitad (para darle tiempo al sistema a que se estabilice)
            ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                
            if (eleccion == 'p') then

                ! Cada valor de magnetización M define una altura en la interfase entre +1 y -1, a su vez una fila en la matriz a partir de la cual se distinguen los dominios. 
                ! Llamo a esa fila 'max_fils_neg', por debajo de la cual hay mayor concentración de -1. Según mis cálculos, es igual a n*densidad de espines negativos:

                max_fils_neg = int(0.5*(n-magn/n))  ! magn = n*n/2 hace que max_fils_neg sea n/4

                if ( (mod(pasos_MC, 4) == 0).AND.(pasos_MC>=fin_pasos*3/4) ) then 
                    
                    ! Calculo la magnetización en los dominios superior (magn_sup, poblada de -1's) e inferior (magn_inf, poblada de +1's)
                    ! Aprovecho para calcular E y M
                    
                    E = 0d0

                    magn_sup = 0
                    do k=1, max_fils_neg
                        do l=1, n
                            magn_sup = magn_sup + s(k,l)
                            E = E + f(k,l)
                        end do  
                    end do

                    magn_inf = 0
                    do k=max_fils_neg+1, n 
                        do l=1, n
                            magn_inf = magn_inf + s(k,l)
                            E = E + f(k,l)
                        end do
                    end do

                    contador = contador + 1
                    promedio_magn_sup = promedio_magn_sup + magn_sup
                    promedio_magn_sup_2 = promedio_magn_sup_2 + magn_sup*magn_sup
                    promedio_magn_inf = promedio_magn_inf + magn_inf


                    promedio_E = promedio_E + E         ! FALTA MULTIPLICAR POR -1/2 !!
                    promedio_EE = promedio_EE + E*E     ! FALTA MULITIPLICAR POR 1/4 !! 
                
    
                end if
            end if

        end do
        ! Finalizan todos los pasos montecarlo

        if (eleccion=='f') close(unit_number)

        if (eleccion=='p') then

            promedio_magn_sup = promedio_magn_sup/contador              ! <M_sup>
            promedio_magn_inf = promedio_magn_inf/contador              ! <M_inf>
            promedio_E = -0.5d0*promedio_E/contador                     ! <E>
            promedio_EE = 0.25d0*promedio_EE/contador                   ! <E^2>
            promedio_magn_sup_2 = promedio_magn_sup_2/contador          ! <M_sup^2> 
            ! promedio_magn_inf_2 = promedio_magn_inf_2/contador          ! <M_inf^2>       ji_inf = ji_sup para M=0

            cn = (promedio_EE-promedio_E*promedio_E)/(n*n*T(indice_T)*T(indice_T))
            ji = (promedio_magn_sup_2-promedio_magn_sup*promedio_magn_sup)/(n*n*T(indice_T))

            ! Encuentro la temperatura crítica en ambas gráficas como el máximo de cn y ji
            if ( cn>=max_cn ) then
                max_cn = cn
                Tc_cn = T(indice_T)
            end if

            if ( ji>=max_ji ) then
                max_ji = ji
                Tc_ji = T(indice_T)
            end if

            write(unit_number2, *) T(indice_T), promedio_magn_sup/(n*n), promedio_magn_inf/(n*n), promedio_E/(n*n), &
            cn, ji
            

            close(unit_number2)
        end if


    end subroutine
    

    ! f(n,m) = s(n,m)*( s(n,m+1) + s(n,m-1) + s(n+1,m) + s(n-1,m) )

    integer*8 function f(k, l)
        integer :: k, k_up, k_down, l, l_right, l_left

        k_up = k+1
        k_down = k-1
        l_right = l+1
        l_left = l-1       

        ! if (k==n) k_up = 1
        ! if (k==1) k_down = n
        if (l==n) l_right = 1
        if (l==1) l_left = n

        ! Esta expresión se utiliza para \Delta E, y excluye a las filas 1 y n pues estas no pueden seleccionarse al no poder cambiar
        f = s(k,l)*( s(k, l_right) + s(k, l_left) + s(k_up, l) + s(k_down, l) )

        ! Estas dos expresiones se usan para incluir las filas 1 y n en el cálculo de E(C)
        if ( k==1 ) f = s(k,l)*( s(k, l_right) + s(k, l_left) + s(k_up, l) )
        if ( k==n ) f = s(k,l)*( s(k, l_right) + s(k, l_left) + s(k_down, l) )

    end function

    ! Función que calcula la diferencia de energía entre dos configuraciones (H(C')-H(C)=delta_E) en el modelo de Kawasaki. Está en función de la expresión
    ! que se utiliza en el cálculo de delta_E para la dinámica de Glauber.

    integer*8 function g(k, l, k_prima, l_prima)
        integer :: k, l, k_prima, l_prima
        g = 2*f(k,l) + 2*f(k_prima, l_prima) - 4*s(k,l)*s(k_prima,l_prima)
    end function

end program kawasaki