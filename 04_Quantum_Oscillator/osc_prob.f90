program osc
    implicit none

    integer :: S, L, n, i, j, ene
    real*8 :: delta_x, delta_t, pi, normalizacion, omega, suma, sqrt_norm
    double complex :: conjugado
    double complex, allocatable :: fi(:), a_2(:), alfa(:), beta(:), b(:), q(:)
    real*8, allocatable :: V(:), V_virgulilla(:)

    real*8 :: tiempo_max, frec, gamma, x_0, factorial, pc

    ! ******************
    ! PROGRAMA PRINCIPAL 
    ! ******************

    ! Parámetros iniciales

    S = 5000
    ! LA FUNCION DE ONDA ES UN VECTOR COMPLEJO DE DIMENSIÓN 1000 QUE VA EVOLUCIONANDO (CADA ELEMENTO ES UN NUM COMPLEJO (dim 2000), igual al valor fi(j*delta_x) para cada j)
    allocate( fi(0:S), V(0:S), V_virgulilla(0:S), a_2(0:S-1), alfa(0:S-1), beta(0:S-1), b(0:S), q(0:S) )
    
    L = 1
    delta_x = dble(L)/dble(S)
    delta_t = 1.d-3
    tiempo_max = 0.1
    frec = 200d0
    gamma = sqrt((frec/2.d0))                   ! NO CONFUNDIR CON ALFA ; en la notación del Bransden es un alfa, pero le cambio el nombre a gamma
    x_0 = 0.5d0                                 ! El punto donde se centra el potencial armónico

    omega = delta_t/delta_x/delta_x
    alfa(S-1) = complex(0.d0, 0.d0)
    beta(S-1) = complex(0.d0, 0.d0)
    q(0) = complex(0.d0, 0.d0)
    pi = acos(-1.d0)

    ! PARA CAMBIAR DE AUTOFUNCIÓN (ene = 0, 3, 10)
    ene = 0
    ! PARA CAMBIAR DE AUTOFUNCIÓN (ene = 0, 3, 10)

    factorial = 1.d0
    i = ene
    do while(i>1)
        factorial = factorial*i
        i = i-1
    end do
    
    normalizacion = sqrt(gamma/(sqrt(pi)*2.d0**ene*factorial ))

    call init_fi(ene) ! Input: ene -> es el parámetro que controla qué autofunción se grafica (n=1,2,3,...)
    call init_V()
    call calculo_de_alfa()

    open(unit=13, file='prob_n0.dat', status='unknown')
    open(unit=13+3, file='prob_n3.dat', status='unknown')
    open(unit=13+25, file='prob_n25.dat', status='unknown')

    do n=0, int(tiempo_max/delta_t)

        ! Calculo la densidad de probabilidad en x y t, la parte real e imaginaria de la función de onda y compruebo la normalización

        suma = 0.d0
        i=0
        do j=0, S
            sqrt_norm = fi(j)*conjg(fi(j))
            suma = suma + sqrt_norm*delta_x
            
            !pc = 1.d0/(pi*sqrt( 2d0*ene+1d0 - (gamma*(j*delta_x-x_0))**2 ))
            pc = 1.d0/sqrt(4.d0*(ene+0.5d0)/frec-(j*delta_x-x_0)**2)/pi
            
            write(13+ene,*) j*delta_x, sqrt_norm, pc
        end do
    

        write(13+ene,*)
        write(13+ene,*)

        ! Evolución: se calcula la función compleja en un tiempo posterior

        call calculo_de_beta()
        call evoluciona_q()
        call evoluciona_fi()

    end do

    close(13)

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
        else if (estado_n==10) then
            hermite = +1024d0*x**10 - 23040d0*x**8 &
            +161280d0*x**6 - 403200d0*x**4 &
            +302400d0*x**2 - 30240d0*1
        else if (estado_n==25) then
            hermite = 33554432d0*x**25 - 5033164800d0*x**23 + 318347673600d0*x**21 - 11142168576000d0*x**19 + &
            238163853312000d0*x**17 - 3239028405043200d0*x**15 + 28341498544128000d0*x**13 - 157902634745856000d0*x**11 &
            + 542790306938880000d0*x**9 - 1085580613877760000d0*x**7 + 1139859644571648000d0*x**5 - 518118020259840000d0*x**3 &
            + 64764752532480000d0*x            
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