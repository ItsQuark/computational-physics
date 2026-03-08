program osc
    implicit none

    integer :: S, L, n, i, j, ene
    real*8 :: delta_x, delta_t, pi, normalizacion, omega, suma, sqrt_norm, suma_x, valormedioXcuadrado, inc_x, inc_p
    double complex :: suma_p, suma_h, conjugado
    double complex, allocatable :: fi(:), a_2(:), alfa(:), beta(:), b(:), q(:)
    real*8, allocatable :: V(:), V_virgulilla(:)

    real*8 :: tiempo_max, frec, gamma, x_0, factorial, media_H

    ! ******************
    ! PROGRAMA PRINCIPAL 
    ! ******************

    ! Parámetros iniciales

    S = 1000
    ! LA FUNCION DE ONDA ES UN VECTOR COMPLEJO DE DIMENSIÓN 1000 QUE VA EVOLUCIONANDO (CADA ELEMENTO ES UN NUM COMPLEJO (dim 2000), igual al valor fi(j*delta_x) para cada j)
    allocate( fi(0:S), V(0:S), V_virgulilla(0:S), a_2(0:S-1), alfa(0:S-1), beta(0:S-1), b(0:S), q(0:S) )
    
    L = 1
    delta_x = dble(L)/dble(S)
    delta_t = 1.d-3
    tiempo_max = 0.5
    frec = 200d0
    gamma = sqrt((frec/2.d0))                   ! NO CONFUNDIR CON ALFA ; en la notación del Bransden es un alfa, pero le cambio el nombre a gamma
    x_0 = 0.5d0                                 ! El punto donde se centra el potencial armónico

    omega = delta_t/delta_x/delta_x
    alfa(S-1) = complex(0.d0, 0.d0)
    beta(S-1) = complex(0.d0, 0.d0)
    q(0) = complex(0.d0, 0.d0)
    pi = acos(-1.d0)

    ! PARA CAMBIAR DE AUTOFUNCIÓN (ene = 0, 1, 2, ...)
    ene = 0

    factorial = 1.d0
    i = ene
    do while(i>1)
        factorial = factorial*i
        i = i-1
    end do
    
    normalizacion = sqrt(gamma/(sqrt(pi)*2.d0**ene*factorial))

    call init_fi(ene) ! Input: ene -> es el parámetro que controla qué autofunción se grafica
    call init_V()
    call calculo_de_alfa()

    ! a_2(j) es constante para todo j ; alfa(j) no hay que calcularlo más veces en todo el programa ; ambos están definidos desde j=0 hasta S-1
    ! Tanto b(j) como beta(j) se va actualizando

    open(unit=10, file='dens_prob.dat', status='unknown')
    open(unit=11, file='normalizado.dat', status='unknown')          
    open(unit=12, file='valores_medios.dat', status='unknown')

    do n=0, int(tiempo_max/delta_t)

        ! Calculo la densidad de probabilidad en x y t, la parte real e imaginaria de la función de onda y compruebo la normalización

        suma = 0.d0
        do j=0,S
            sqrt_norm = fi(j)*conjg(fi(j))
            suma = suma + sqrt_norm*delta_x
            write(10,*) j*delta_x, sqrt_norm, real(fi(j)), aimag(fi(j))
        end do

        write(10,*)
        write(10,*)

        write(11,*) n, suma
        write(11,*)
        write(11,*)

        ! Calculo los valores medios de x, p y H y comparamos con el valor teórico   
        
        ! <x> = suma_x      <p> = suma_p        <H_caja> = suma_h = <p^2> -> media_H = <H_osc> = <p^2> + frec^2/4 * ( <x^2> - 2*x_0*<x> + x_0^2 )
        ! <x^2> = valormedioXcuadrado       <p^2> = suma_h

        suma_x = 0.d0
        valormedioXcuadrado = 0.d0
        suma_p = complex(0.d0, 0.d0)
        suma_h = complex(0.d0, 0.d0)

        do j=0,S
            conjugado = conjg(fi(j))
            suma_x = suma_x + j*fi(j)*conjugado
            if (j.ne.S) suma_p = suma_p + conjugado*(fi(j+1)-fi(j))
            if ((j.ne.0).or.(j.ne.S)) suma_h = suma_h + conjugado*( fi(j+1) - 2.d0*fi(j) + fi(j-1) )
            valormedioXcuadrado = valormedioXcuadrado + j*j*fi(j)*conjugado
        end do

        suma_x = delta_x*delta_x*suma_x
        suma_p = complex(0.d0, -1.d0)*suma_p
        suma_h = -1.d0/delta_x*suma_h

        valormedioXcuadrado = delta_x*delta_x*delta_x*valormedioXcuadrado   ! Aquí finaliza el cálculo de <x^2>, por lo que después calculamos <H>

        ! Izq: valor teórico    Der: programa                        Todo coincide :)

        !print*, x_0, suma_x                                                
        !print*, (2*ene+1.d0)/frec + x_0*x_0, valormedioXcuadrado
        !print*, 0.d0, real(suma_p)           
        !print*, (ene+1.d0/2.d0)*frec/2.d0, real(suma_h)   
        !print*, ene+0.5d0, inc_p*inc_x,                                                       

        media_H = real(suma_h) + frec*frec/4.d0*(valormedioXcuadrado-x_0*x_0)

        inc_x = sqrt(valormedioXcuadrado-suma_x*suma_x)
        inc_p = sqrt(real(suma_h)-real(suma_p)*real(suma_p))

        write(12,*) n, suma_x, real(suma_p), media_H, (ene+0.5d0)*frec, inc_x*inc_p, ene + 0.5d0
        write(12,*)


        ! Evolución: se calcula la función compleja en un tiempo posterior

        call calculo_de_beta()
        call evoluciona_q()
        call evoluciona_fi()

    end do

    close(10)
    close(11)
    close(12)

    stop
    contains

    ! **********
    ! SUBRUTINAS
    ! **********

    ! Función: genera los polinomios de Hermite

    real*8 function hermite(estado_n, x)
        real*8 :: x
        integer :: estado_n

        if (estado_n==0) then
            hermite = 1.d0
        else if (estado_n==1) then
            hermite = 2.d0*x
        else if (estado_n==2) then
            hermite = 4.d0*x*x-2.d0
        else if (estado_n==3) then
            hermite = 8d0*x*x*x - 12d0*x
        end if
    end function


    ! Función: genera los S valores iniciales de la autofunción del hamiltoniano en el vector fi

    real*8 function f(x,estado_n)
        real*8 :: x
        integer :: estado_n
        f = normalizacion*exp(-gamma*gamma*(x-x_0)**2/2.d0)*hermite(estado_n, gamma*(x-x_0))
    end function

    ! Subrutina: inicializa la función de onda (vector complejo de dim 1000). Recibe la autofunción n-ésima (n = 1,2,3) correspondientemente normalizada
    !            aprovecho para inicializar el vector b (en el tiempo inicial), necesario para calcular beta

    subroutine init_fi(estado_n)
            implicit none
            integer :: j, estado_n
            do j = 0, S
                fi(j) = complex(f(j*delta_x, estado_n), 0.d0)
            end do
    end subroutine init_fi

    ! Subrutina: inicializa el potencial V(x), el potencial virgulilla, el vector a_2,j (usado para calcular alfa)

    subroutine init_V()
        implicit none
        integer :: j
        do j = 0, S
            V(j) = frec*frec/4.d0*(j*delta_x-x_0)**2
        end do
        V_virgulilla = delta_x*delta_x*V

        do j=0,S-1
            a_2(j) = complex(-2.d0-V_virgulilla(j), 2.d0/omega)
        end do

    end subroutine init_V

    ! Subrutina: se calcula alfa (vector constante en todo el algoritmo porque V no cambia en el tiempo)

    subroutine calculo_de_alfa()
        implicit none
        integer :: j
        do j=S-1, 0, -1
            alfa(j-1) = -1.d0/(a_2(j)+alfa(j))
        end do
    end subroutine

    ! Subrutina: se calcula beta en un tiempo concreto n
    
    subroutine calculo_de_beta()
        implicit none
        integer :: j
        do j = 0, S
            b(j) = complex(0.d0, 4.d0/omega)*fi(j)
        end do

        do j = S-1, 0, -1
            beta(j-1) = (b(j) - beta(j) )/( a_2(j) + alfa(j) )
        end do
    end subroutine

    ! Subrutina: se calcula q(j) para un tiempo n

    subroutine evoluciona_q()
        implicit none
        integer :: j

        do j = 0, S-1
            q(j+1) = alfa(j)*q(j) + beta(j)
        end do

    end subroutine

    ! Subrutina: se calcula phi en un tiempo posterior

    subroutine evoluciona_fi()
        implicit none
        integer :: j

        do j = 0, S
            fi(j) = q(j) - fi(j)
        end do
    end subroutine
    
end program osc