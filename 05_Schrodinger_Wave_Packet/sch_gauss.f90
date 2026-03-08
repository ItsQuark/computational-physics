program sch
    implicit none

    integer :: S, L, n, i, j
    real*8 :: delta_x, delta_t, pi, normalizacion, omega, suma, sqrt_norm, suma_x, valormedioXcuadrado, inc_x, inc_p
    double complex :: suma_p, suma_h, conjugado, complex_i
    double complex, allocatable :: fi(:), a_2(:), alfa(:), beta(:), b(:), q(:)
    real*8, allocatable :: V(:), V_virgulilla(:)

    ! ******************
    ! PROGRAMA PRINCIPAL 
    ! ******************

    ! Parámetros iniciales

    S = 1000
    ! LA FUNCION DE ONDA ES UN VECTOR COMPLEJO DE DIMENSIÓN 1000 QUE VA EVOLUCIONANDO (CADA ELEMENTO ES UN NUM COMPLEJO (dim 2000), igual al valor fi(j*delta_x) para cada j)
    allocate( fi(0:S), V(0:S), V_virgulilla(0:S), a_2(0:S-1), alfa(0:S-1), beta(0:S-1), b(0:S), q(0:S) )
    
    L = 1
    delta_x = dble(L)/dble(S)
    delta_t = 1.d-4
    omega = delta_t/delta_x/delta_x
    alfa(S-1) = complex(0.d0, 0.d0)
    beta(S-1) = complex(0.d0, 0.d0)
    q(0) = complex(0.d0, 0.d0)
    
    pi = acos(-1.d0)
    complex_i = complex(0.d0, 1.d0)

    ! Lo calculo una vez y me evito tener que llamar a la función sqrt cada vez que llame a la función fi
    normalizacion = 1d0/sqrt(0.11077836488271293)
    
    call init_fi() ! Input: ene -> es el parámetro que controla qué autofunción se grafica (n=1,2,3,...)
    call init_V()
    call calculo_de_alfa()
 
    ! a_2(j) es constante para todo j ; alfa(j) no hay que calcularlo más veces en todo el programa ; ambos están definidos desde j=0 hasta S-1
    ! Tanto b(j) como beta(j) se va actualizando

    open(unit=10, file='dens_prob_g.dat', status='unknown')                       ! j, dens_prob, parte real, parte imaginaria
    open(unit=11, file='normalizado_g.dat', status='unknown')                     ! n, integral
    open(unit=12, file='valores_medios_g.dat', status='unknown')

    do n=0, 500

        ! Calculo la densidad de probabilidad en x y t, la parte real e imaginaria de la función de onda y compruebo la normalización

        suma = 0.d0
        do j=0,S
            sqrt_norm = fi(j)*conjg(fi(j))
            suma = suma + sqrt_norm*delta_x
            write(10,*) j*delta_x, sqrt_norm, real(fi(j)), aimag(fi(j))
        end do

        print*, suma

        write(10,*)
        write(10,*)

        write(11,*) n, suma
        write(11,*)
        write(11,*)

        ! Calculo los valores medios de x, p y H y comparamos con el valor teórico        

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

        inc_x = sqrt(valormedioXcuadrado-suma_x*suma_x)
        inc_p = sqrt(real(suma_h)-real(suma_p)*real(suma_p))

        write(12,*) n, suma_x, real(suma_p), real(suma_h), inc_x*inc_p
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

    ! Subrutina: inicializa la función de onda (vector complejo de dim 1000). Recibe la autofunción n-ésima (n = 1,2,3) correspondientemente normalizada
    !            aprovecho para inicializar el vector b (en el tiempo inicial), necesario para calcular beta

    subroutine init_fi()
            implicit none
            integer :: j
            do j = 0, S
                ! fi(j) = complex(f(j*delta_x, n), 0.d0)
                fi(j) = normalizacion*exp( -8.d0*(4.d0*j-S)*(4.d0*j-S)/(S*S) )*exp(complex_i/(2.d0*sqrt(delta_t))*j*delta_x)

            end do
    end subroutine init_fi

    ! Subrutina: inicializa el potencial V(x), el potencial virgulilla, el vector a_2,j (usado para calcular alfa)

    subroutine init_V()
        implicit none
        integer :: j
        do j = 0, S
            V(j) = 0.d0
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
    
end program sch