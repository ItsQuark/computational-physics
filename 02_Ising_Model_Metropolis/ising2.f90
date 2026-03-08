program ising
    use randomnumber
    implicit none
    integer :: n
    integer*8 :: resultado, fin_pasos, archivo_append
    integer, allocatable :: s(:,:)
    character(len=1) :: eleccion

    call dran_ini(18627)

    ! Longitud de la red (#puntos = n²)
    n = 128
    allocate (s(n,n))

    
    ! Algoritmo metropolis (input: temperatura) -> Genera los .dat necesarios para hacer el gif o los datos para el promedio.
    
    print*, "¿Desea escribir los datos en un fichero necesario para crear los gifs (f) o calcular los promedios (p)?"
    ! f: por defecto, el algoritmo usa 100 pasos montecarlo debido al peso de los archivos (20 mb) 
    fin_pasos = 100
    ! p: si elige promediar, no se escribe nada en un fichero para no saturarlos, y doy la posibilidad de elegir muchos pasos montecarlo
    read(*,*) eleccion
    if (eleccion == "p") then
        print*, "Pasos montecarlo (numero entero par mayor que cero):"
        read(*,*) fin_pasos
    end if

    ! 1: todos 1 // -1: todos -1 // 0: ordenado, mismos 1 que -1 (salvo n impar, en tal caso hay un +1 neto) // 10: aleatorio
    call init(0)

    archivo_append = 0
    call metropolis(5.d0)
    archivo_append = 1
    call init(0)
    call metropolis(2.27d0)
    call init(0)
    call metropolis(1.d0)
    call init(0)
    call metropolis(0.01d0)

    stop
    contains

    ! Subrutina init: inicializa la red del modelo de Ising a la Configuracion deseada. 1: todos 1 // -1: todos 1 // 0: ordenado, mismos 1 que -1 (salvo n impar, en tal caso hay un +1
    !                 neto) // 10: aleatorio

    subroutine init(estado)
        implicit none
        integer :: estado, k, l

        if (estado == 1) then
            s = 1

        else if (estado == -1) then
            s = -1

        else if (estado == 0) then

            do k=1,n
                do l=1,n
                    if (mod(l+k,2)==0) s(k,l) = 1
                    if (mod(l+k,2)/=0) s(k,l) = -1
                end do
            end do

        else if (estado == 10) then

            ! Implementar aleatorio

        else

        print*, "Introduzca un valor adecuado (1, -1, 0, 10)"

        end if

    end subroutine

    subroutine metropolis(T)
        implicit none
        real*8 :: T, p, promedio, E, promedio_M
        integer*8 :: k, l, q, r, delta_E, iter, pasos_MC, contador, magnetizacion

        ! Ficheros donde almaceno la evolucion de la red de spines para distintas temperaturas, y la energía

        open(unit=10, file="5.dat", status="unknown")
        open(unit=11, file="227.dat", status="unknown")
        open(unit=12, file="1.dat", status="unknown")
        open(unit=13, file="01.dat", status="unknown")
        if (archivo_append == 1) then
            open(unit=14, file="energia.dat", status="unknown", position="append")
        else if (archivo_append == 0) then
            open(unit=14, file="energia.dat", status="unknown")
        end if

        ! Datos iniciales, antes de meterse en el bucle que itera sobre cada paso montecarlo

        promedio = 0.d0
        promedio_M = 0.d0
        contador = 0

        do pasos_MC = 1, fin_pasos


        ! Escribo la evolucion de la matriz para la temperatura deseada en el fichero correspondiente, si el usuario eligió 'f'

        if (eleccion == 'f') then

            if (T==5.d0) then
                do k=1,n
                    write(10,*) s(k,:)
                end do
            write(10,*)
            write(10,*)

        else if (T==2.27d0) then

            do k=1,n
                write(11,*) s(k,:)
            end do
            write(11,*)
            write(11,*)

        else if (T==1.d0) then

            do k=1,n
                write(12,*) s(k,:)
            end do
            write(12,*)
            write(12,*)

        else if (T==0.01d0) then

            do k=1,n
                write(13,*) s(k,:)
            end do
            write(13,*)
            write(13,*)
        end if

        end if

        ! Procedo con un paso montecarlo: se hacen n*n preguntas en casillas aleatorias

        do iter = 1, n*n

        ! Elijo un punto al azar de la red

        k = i_dran(int(n, 8))
        l = i_dran(int(n, 8))

        ! Calculo delta_E

        call termino(k,l)

        delta_E = 2*resultado

        ! Calculo p

        p = min(1.d0, exp(-delta_E/T))

        ! Decisión

        if ( dran_u()<p ) s(k,l) = -s(k,l)

        end do

        ! Finaliza un paso montecarlo, aprovecho este espacio para realizar los cálculos pertinentes a este paso MC, antes de iniciar otro.

        ! Calculo la energía y la magnetización

        if ( (mod(pasos_MC, 5) == 0).AND.(pasos_MC>fin_pasos/2) ) then 

        ! La energía y la magnetización se calcula mediante un sumatorio, lo inicializo a cero
        E = 0.d0
        magnetizacion = 0

        do q=1, n
                do r=1, n
                    call termino(q,r)
                    E = E + resultado
                    magnetizacion = magnetizacion + s(q,r)
                end do
        end do

        E = -0.5*E

        promedio_M = promedio_M + magnetizacion
        promedio = promedio + E
        contador = contador + 1

        end if

        end do

    ! Finalizan todos los pasos montecarlo

    promedio = promedio/(n*n*contador)
    promedio_M = abs(promedio_M)/(n*n*contador)
    write(14,*) T, promedio, promedio_M

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)

    end subroutine
    

    ! Calcula el término que está dentro del sumatorio

    subroutine termino(k, l)
        implicit none
        integer*8 k, l

        if ( (k>1).AND.(k<n).AND.(l>1).AND.(l<n) ) then
            resultado = s(k,l)*( s(k+1,l) + s(k-1,l) + s(k, l+1) + s(k, l-1)  )

        else if ( (k==1).AND.(l==1) ) then
            resultado = s(k,l)*( s(k+1,l) + s(n,l) + s(k, l+1) + s(k, n)  )

        else if ( (k==1).AND.(l>1).AND.(l<n) ) then
            resultado = s(k,l)*( s(k+1,l) + s(n,l) + s(k, l+1) + s(k, l-1)  )

        else if ( (k==1).AND.(l==n) ) then
            resultado = s(k,l)*( s(k+1,l) + s(n,l) + s(k, 1) + s(k, l-1)  )

        else if ( (l==n).AND.(k>1).AND.(k<n) ) then
            resultado = s(k,l)*( s(k+1,l) + s(k-1,l) + s(k, 1) + s(k, l-1)  )

        else if ( (k==n).AND.(l==n) ) then
            resultado = s(k,l)*( s(1,l) + s(k-1,l) + s(k, 1) + s(k, l-1)  )

        else if ( (k==n).AND.(l>1).AND.(l<n) ) then
            resultado = s(k,l)*( s(1,l) + s(k-1,l) + s(k, l+1) + s(k, l-1)  )

        else if ( (k==n).AND.(l==1) ) then
            resultado = s(k,l)*( s(1,l) + s(k-1,l) + s(k, l+1) + s(k, n)  )

        else if ( (l==1).AND.(k>1).AND.(k<n) ) then
            resultado = s(k,l)*( s(k+1,l) + s(k-1,l) + s(k, l+1) + s(k, n)  )

        end if

    end subroutine

end program
