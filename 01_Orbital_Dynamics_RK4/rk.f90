program rk
    implicit none

    integer :: i, contador, lineas
    real*8 :: constante_G, m_t, m_l, d_tl, omega, r_t, r_l, delta, mu, t, tiempo_maximo, pi, v_orbital, h, phi_inicial, dias, v
    real*8 :: theta_inicial, m, hamiltoniano
    real*8, allocatable :: y(:), k1(:), k2(:), k3(:), k4(:)

    ! Datos del problema

    constante_G = 6.67d-11
    m_t = 5.9736d24
    m_l = 0.07349d24
    d_tl = 3.844d8
    omega = 2.6617d-6
    r_t = 6.378160d6
    r_l = 1.7374d6
    pi = acos(-1.d0)

    m = 1.d0    ! ???

    delta = constante_G*m_t/(d_tl*d_tl*d_tl)
    mu = m_l/m_t

    ! Construyo los vectores que almacenan las coordenadas generalizadas y las distintas k's. OJO: ki = k_super_i, no confundir con índice de coordenada generalizada.
    ! El orden es importante: 1 -> r    2 -> phi    3 -> p_r    4 -> p_phi

    allocate(y(1:4), k1(1:4), k2(1:4), k3(1:4), k4(1:4))

    t = 0.d0
    
    ! Zona de control

    h = 50
    dias = 8
    tiempo_maximo = dias*24*3600

    phi_inicial = pi/2.d0
    theta_inicial = 0.d0                            ! Interactúa en 0.95 rad

    v_orbital = sqrt(constante_G*m_t/r_t)
    v = v_orbital*sqrt(2.d0)

    contador = 0
    lineas = 0

    call init(phi_inicial, theta_inicial, v_orbital)

    print*, y

    open(unit = 11, file = "datos.dat", status = "unknown")
    open(unit = 12, file = "constante.dat", status = "unknown")
    
    do while (t<=tiempo_maximo)

        if (mod(contador, 2)==0) then
            lineas = lineas + 1
            write(11, *) cos(omega*t), sin(omega*t), y(1)*cos(y(2)), y(1)*sin(y(2)), t/3600.d0
            write(11, *)
            write(11, *)
        end if 
        
        ! Constante del movimiento

        !hamiltoniano = -omega*m*d_tl*d_tl*y(4) + 0.5d0*m*d_tl*d_tl*y(3)*y(3) + 0.5d0*m*d_tl*d_tl*(y(4)/y(1))**2 - constante_G*m*m_t/(y(1)*d_tl)

        !write(12,*) t/3600.d0, hamiltoniano - constante_G*m*m_l/sqrt(d_tl*d_tl + y(1)*y(1)*d_tl*d_tl - 2.d0*y(1)*d_tl*d_tl*cos(y(2)-omega*t))

        hamiltoniano = -omega*m*d_tl*d_tl*y(4) + 0.5d0*m*d_tl*d_tl*y(3)*y(3) + &
               0.5d0*m*d_tl*d_tl*(y(4)/y(1))**2 - constante_G*m*m_t/(y(1)*d_tl)
        
        write(12,*) t/3600.d0, hamiltoniano - constante_G*m*m_l/ &
             sqrt(d_tl*d_tl + y(1)*y(1)*d_tl*d_tl - 2.d0*y(1)*d_tl*d_tl*cos(y(2)-omega*t))


        ! Evaluamos las k's

        do i = 1, 4
            k1(i) = h*f(i, y(1), y(2), y(3), y(4), t)
        end do

        do i = 1, 4
            k2(i) = h*f(i, y(1)+k1(1)/2.d0, y(2)+k1(2)/2.d0, y(3)+k1(3)/2.d0, y(4)+k1(4)/2.d0, t+h/2.d0)
        end do

        do i = 1, 4
            k3(i) = h*f(i, y(1)+k2(1)/2.d0, y(2)+k2(2)/2.d0, y(3)+k2(3)/2.d0, y(4)+k2(4)/2.d0, t+h/2.d0)
        end do

        do i = 1, 4
            k4(i) = h*f(i, y(1)+k3(1), y(2)+k3(2), y(3)+k3(3), y(4)+k3(4), t+h)
        end do

        y = y + 1.d0/6.d0*(k1 + 2*k2 + 2*k3 + k4)
        t = t + h
        contador = contador + 1

    end do

    print*, contador, "iteraciones"
    print*, lineas, "líneas escritas"
    print*, v/1000.d0, "km/s"
    
    close(11)
    close(12)

    contains

    real*8 function f(n, y1, y2, y3, y4, tiempo)
        integer :: n
        real*8 :: y1, y2, y3, y4, tiempo, y1_prima
        if (n==1) then
            f = y3
        else if (n==2) then
            f = y4/(y1*y1)
        else if (n==3) then
            y1_prima = sqrt( 1.d0 + y1*y1 - 2.d0*y1*cos(y2-omega*tiempo) )
            f = y4*y4/(y1*y1*y1) - delta*( 1.d0/(y1*y1) + mu/(y1_prima*y1_prima*y1_prima)*(y1-cos(y2-omega*tiempo)) )
        else if (n==4) then
            y1_prima = sqrt( 1.d0 + y1*y1 - 2.d0*y1*cos(y2-omega*tiempo) )
            f = -delta*mu*y1/(y1_prima*y1_prima*y1_prima)*sin(y2-omega*tiempo)
        end if
    end function

    ! Subrutina: inicializa el vector de coordenadas y momentos generalizados 'y' ; recibe de input las CI SIN reescalar

    subroutine init(phi_0, theta_0, v_0)    
        implicit none
        ! r_0 siempre estará localizado en la superficie terrestre; el resto son de libre configuración: phi_0 es la latitud, theta_0 el ángulo de lanzamiento
        ! del cohete en la superficie con respecto del eje x, v_0 es el módulo de la velocidad de lanzamiento (aprox. la vel. de escape)
        real*8 :: phi_0, theta_0, v_0

        y = [r_t/d_tl, phi_0, v_0/d_tl*cos(theta_0-phi_0), r_t/d_tl*v_0/d_tl*sin(theta_0-phi_0)]    ! Se introducen las variables y se aplica el reescalamiento
    end subroutine init

end program rk