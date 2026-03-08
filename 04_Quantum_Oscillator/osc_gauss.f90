program osc_gauss
    implicit none

    integer :: S, L, n, i, j, ene
    real*8 :: delta_x, delta_t, pi, normalizacion, omega, suma, sqrt_norm, suma_x, valormedioXcuadrado, inc_x, inc_p
    double complex :: suma_p, suma_h, conjugado, complex_i
    double complex, allocatable :: fi(:), a_2(:), alfa(:), beta(:), b(:), q(:)
    real*8, allocatable :: V(:), V_virgulilla(:)

    real*8 :: tiempo_max, frec, gamma, x_0, media_H, pos_clasica, E_clasico, c_0, mom_clasico, xbar_0, pbar_0, desv_est

    ! ******************
    ! PROGRAMA PRINCIPAL 
    ! ******************

    ! Parámetros iniciales

    S = 1000
    ! LA FUNCION DE ONDA ES UN VECTOR COMPLEJO DE DIMENSIÓN 1000 QUE VA EVOLUCIONANDO (CADA ELEMENTO ES UN NUM COMPLEJO (dim 2000), igual al valor fi(j*delta_x) para cada j)
    allocate( fi(0:S), V(0:S), V_virgulilla(0:S), a_2(0:S-1), alfa(0:S-1), beta(0:S-1), b(0:S), q(0:S) )
    
    L = 1.d0
    delta_x = dble(L)/dble(S)
    delta_t = 1.d-3
    tiempo_max = 0.5                            ! Guion pide: 0.5
    frec = 200d0
    gamma = sqrt((frec/2.d0))                   ! NO CONFUNDIR CON ALFA ; en la notación del Bransden es un alfa, pero le cambio el nombre a gamma
    x_0 = 0.5d0                                 ! El punto donde se centra el potencial armónico

    omega = delta_t/delta_x/delta_x
    alfa(S-1) = complex(0.d0, 0.d0)
    beta(S-1) = complex(0.d0, 0.d0)
    q(0) = complex(0.d0, 0.d0)
    pi = acos(-1.d0)
    complex_i = complex(0.d0, 1.d0)

    c_0 = 0.3d0
    desv_est = 1d0/10d0
    normalizacion = 1.d0/sqrt(0.17724538509029283)!1.d0/sqrt(0.11077836568101400) !
    print*, normalizacion
    

    call init_fi(c_0, desv_est)
    call init_V()
    call calculo_de_alfa()
 
    ! a_2(j) es constante para todo j ; alfa(j) no hay que calcularlo más veces en todo el programa ; ambos están definidos desde j=0 hasta S-1
    ! Tanto b(j) como beta(j) se va actualizando

    open(unit=10, file='dens_prob_g.dat', status='unknown')
    open(unit=11, file='normalizado_g.dat', status='unknown')          
    open(unit=12, file='valores_medios_g.dat', status='unknown')

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

        valormedioXcuadrado = delta_x*delta_x*delta_x*valormedioXcuadrado   

        ! Izq: valor teórico    Der: programa       esto para t=0                   Todas coinciden! Por tanto, <H> y la incertidumbre están bien

        !print*, <x> = 0.2999999644813863, suma_x = 0.29999997224902647                                              
        !print*, <x^2> = 0.09195311, valormedioXcuadrado = 0.09195311
        !print*, <p> = 0, real(suma_p) = 0           
        !print*, <p^2> = 127.9999881929782, real(suma_h) = 127.99589310246206       OJO: real(suma_h) = <p^2>                                                

        media_H = real(suma_h) + frec*frec/4.d0*(valormedioXcuadrado+x_0*x_0-2d0*suma_x*x_0)

        ! Calculamos dos valores que determinan las condiciones iniciales de la ecuación de movimiento para el oscilador clásico, asegurando que E_q = E_clasico

        if (n==0) then
            xbar_0 = sqrt(valormedioXcuadrado-2d0*x_0*suma_x+x_0*x_0)
            pbar_0 = sqrt(real(suma_h))
        end if

        pos_clasica = x_0 + xbar_0*cos(frec*n*delta_t) + 2d0/frec*pbar_0*sin(frec*n*delta_t)
        mom_clasico = -0.5d0*frec*xbar_0*sin(frec*n*delta_t) + pbar_0*cos(frec*n*delta_t)
        E_clasico = mom_clasico*mom_clasico + frec*frec/4.d0*(pos_clasica-x_0)**2

        if (n==0) then
            print*, "x_clasica, p_clasico, E_clasico", pos_clasica, mom_clasico, E_clasico
            print*, suma_x, real(suma_p), media_H
        end if

        inc_x = sqrt(valormedioXcuadrado-suma_x*suma_x)
        inc_p = sqrt(real(suma_h)-real(suma_p)*real(suma_p))

        write(12,*) n, suma_x, real(suma_p), media_H, inc_x*inc_p, pos_clasica, mom_clasico, E_clasico
        write(12,*)

        if (n==0) then
            print*, "<x>=", suma_x, "<x^2>=", valormedioXcuadrado
            print*, "<p>=", real(suma_p), "<p^2>=", real(suma_h)
            print*, "<H>=", media_H, "Δx·Δp=", inc_x*inc_p
        end if


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

    ! Subrutina: inicializa la función de onda (vector complejo de dim 1000). Recibe la autofunción n-ésima (n = 1,2,3) correspondientemente normalizada
    !            aprovecho para inicializar el vector b (en el tiempo inicial), necesario para calcular beta

    subroutine init_fi(c_0, desv_tip)
        implicit none
        integer :: j
        real*8 :: c_0, desv_tip
        do j = 0, S
            fi(j) = normalizacion*exp( -((j*delta_x-c_0)**2)/(2d0*desv_tip**2) )
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
    
end program osc_gauss