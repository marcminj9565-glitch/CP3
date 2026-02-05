program ising2d
    implicit none

    integer, parameter :: L = 20
    integer, parameter :: N = L*L
    integer :: i, j, step, mc_steps
    integer :: ip, im, jp, jm
    integer :: dE
    real :: T, beta, r
    integer :: spin(L, L)
    real :: E, M
    real :: J

    call random_seed()

    J = 1.0
    T = 2.3
    beta = 1.0 / T
    mc_steps = 100000

    ! Initialize spins randomly
    do i = 1, L
        do j = 1, L
            call random_number(r)
            if (r < 0.5) then
                spin(i,j) = -1
            else
                spin(i,j) = 1
            end if
        end do
    end do

    ! Monte Carlo simulation
    do step = 1, mc_steps
        do i = 1, L
            do j = 1, L

                ip = mod(i, L) + 1
                im = mod(i-2, L) + 1
                jp = mod(j, L) + 1
                jm = mod(j-2, L) + 1

                dE = 2 * J * spin(i,j) * &
                     ( spin(ip,j) + spin(im,j) + &
                       spin(i,jp) + spin(i,jm) )

                if (dE <= 0) then
                    spin(i,j) = -spin(i,j)
                else
                    call random_number(r)
                    if (r < exp(-beta * dE)) then
                        spin(i,j) = -spin(i,j)
                    end if
                end if

            end do
        end do
    end do

    ! Calculate energy and magnetization
    E = 0.0
    M = 0.0

    do i = 1, L
        do j = 1, L
            ip = mod(i, L) + 1
            jp = mod(j, L) + 1

            E = E - J * spin(i,j) * (spin(ip,j) + spin(i,jp))
            M = M + spin(i,j)
        end do
    end do

    print *, "Temperature =", T
    print *, "Energy per spin =", E / N
    print *, "Magnetization per spin =", M / N

end program ising2d
