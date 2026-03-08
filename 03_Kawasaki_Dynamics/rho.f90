program kawasaki
    use randomnumber
    implicit none
    integer :: n, dim_T, w, v, magn
    integer*8 :: fin_pasos
    integer, allocatable :: s(:,:)
    real*8, allocatable :: T(:)
    character(len=1) :: eleccion

    ! README: este programa es una copia modificada de kawasaki.f90, por lo que hay muchas cosas a ignorar. El objetivo es crear las gráficas del promedio de ρ(y)
    ! para cuatro temperaturas distintas, para N = {32, 64, 128}. Se pueden modificar los siguientes parámetros en el bloque de código que encierra el 'else if (-
    ! eleccion=="p") then' en la línea 31: temperaturas del vector T ; magnetización en la variable 'magn'. Luego todos los datos se representan con 'rho.plt'.

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


    ! Distingo la funcionalidad del ejecutable: generación de gifs (fuerzo pocos pasos MC) o cálculo de promedios (doy la posibilidad de dar más pasos MC)
    
    print*, "Presione (p)"
    ! f: por defecto, el algoritmo usa 100 pasos montecarlo debido al peso de los archivos (100 pasos MC = 18.76 MB) 
    fin_pasos = 100
    ! p: si elige promediar, no se escribe nada en un fichero para no saturarlos, y doy la posibilidad de elegir muchos pasos montecarlo
    read(*,*) eleccion
    if (eleccion == "p") then
        print*, "Pasos montecarlo (numero entero par mayor que cero):"
        read(*,*) fin_pasos
    else
        print*, "Error"
        stop
    end if

    call dran_ini(1968339)

    if ( eleccion=="f" ) then

        ! Si solo se desea hacer los gifs, ejecutamos esta parte del código en la que la longitud de red n permanece fija
        ! Longitud de la red (# puntos = n²) (solo funciona para n par, de modo que se permita una magnetización nula)
        n = 64
        magn = 0

        ! Creación del vector que almacena las temperaturas. Para la elección 'f', son 4 cuatro temperaturas
        dim_T = 4
        allocate (T(dim_T))
        T = [1d0, 2d0, 3d0, 10d0]   

        allocate (s(n,n))
        do w = 1, dim_T

            ! Llamada a la subrutina que inicializa la matriz con la magnetización M deseada. Debe ser par, de otro modo no existe el número de espines positivos tal que de esa magnetización impar
            ! La magnetización puede tomar los valores desde -n*n+2*n a n*n-2*n, siendo par, para asegurar que haya el suficiente número de espines para llenar
            call init(magn)
    

            ! Llamada al algoritmo metropolis (input: índice w que etiqueta la temperatura deseada) tantas veces como temperaturas haya
            call metropolis(w)
        end do
    
    else if (eleccion=="p") then
        ! Si se desea hacer cálculo se ejecuta esta parte, donde varían tanto la temperatura como n=[32, 64, 128] a M/(n*n) = m fijo
        n=32
        magn = n*n/2   ! Si M distinto de 0, decomentar
        !magn = 0        ! Si M = 0, decomentar

        dim_T = 4
        allocate(T(dim_T))
        T = [0.01d0, 1.5d0, 2d0, 8d0]

        do v = 1, 3
            allocate (s(n,n)) 
            do w = 1, dim_T
                call init(magn)
                call metropolis(w)
                print*, "Temperatura", T(w), "hecho"
            end do
            print*, "Longitud de red", n, "hecho"

            n = n+n
            magn = n*n/2
            !magn = 0
            deallocate(s)
        end do
        
    else 
        print*, "Solo puede insertar 'f' o 'p'"
    end if
    
    stop
    contains

    ! Subrutina init: inicializa la red del modelo de Ising con la magnetización M deseada, con spin -1 en la fila 1, y +1 en la fila N
    ! M va desde -N*N a + N*N, debe ser par y mayor o igual que 2*n-n*n para que haya suficientes espines para llenar la fila 1

    subroutine init(magnetizacion)
        implicit none
        integer :: k, l, magnetizacion, espines_positivos, desorden, pos11, pos12, pos21, pos22, aux

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

        ! Mezclo las orientaciones de los espines intercambiando posiciones de forma aleatoria, 'desorden' veces
        ! desorden = 10
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
        integer :: indice_T, unit_number, k, l, vecino, k_prima, l_prima, aux, unit_number2, unit_number3
        character(len=200) :: filename, filename2, filename3

        integer*8 :: suma
        real*8, allocatable :: suma_fila(:)

        allocate(suma_fila(n))
 
        unit_number3 = 30+dim_T+v
        write(filename2, *) indice_T
        write(filename3, *) n
        filename3 = trim(adjustl(filename3)) // '_' // trim(adjustl(filename2)) // 'rho.dat'

        ! Ficheros donde se almacenan para distintas temperatuas la evolución de la red de espines (1.dat, 2.dat, ...)
        if (eleccion=='f') open(unit=unit_number, file=filename, status="unknown")
        
        ! Ficheros donde se almacenan para cada n fija, variando T: rho(y)
        if (eleccion=='p') then
            open(unit=unit_number3, file=filename3, status="unknown")       
        end if

        
        ! Antes de iniciar el bucle que itera sobre pasos MC, se inicializan las variables que almacenan promedios

        contador = 0
        suma_fila = 0
        
        do pasos_MC = 1, fin_pasos

            ! Escribo la evolucion de la matriz para la temperatura deseada en el fichero correspondiente, si el usuario eligió 'f'

            if (eleccion == 'f') then
                do k=n, 1, -1
                    ! La escribo al revés porque por razones que desconozco Gnuplot representa los datos 'del revés'
                    write(unit_number,*) s(k,:)
                end do
                write(unit_number,*)
                write(unit_number,*)
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

                

                if ( (mod(pasos_MC, 4) == 0).AND.(pasos_MC>=fin_pasos*3/4) ) then 
                    
                    ! Voy acumulando en el elemento y-ésimo del vector 'suma_fila' los valores de la suma de los espines de la fila y-ésima para luego hacer la media

                    contador = contador + 1

                    do k = 1, n
                        suma = 0
                        do l=1,n
                            suma = suma + s(k, l)
                        end do
                        suma_fila(k) = suma_fila(k) + suma
                    end do
    
                end if
            end if

        end do
        ! Finalizan todos los pasos montecarlo    

        if (eleccion=='f') close(unit_number)
        if (eleccion=='p') then

            suma_fila = suma_fila/n/contador
            
            do k = 1,n
                write(unit_number3, *) k, suma_fila(k)
            end do

            close(unit_number3)
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

        f = s(k,l)*( s(k, l_right) + s(k, l_left) + s(k_up, l) + s(k_down, l) )

    end function

    ! Función que calcula la diferencia de energía entre dos configuraciones (H(C')-H(C)=delta_E) en el modelo de Kawasaki. Está en función de la expresión
    ! que se utiliza en el cálculo de delta_E para la dinámica de Glauber.

    integer*8 function g(k, l, k_prima, l_prima)
        integer :: k, l, k_prima, l_prima
        g = 2*f(k,l) + 2*f(k_prima, l_prima) - 4*s(k,l)*s(k_prima,l_prima)
    end function

end program kawasaki