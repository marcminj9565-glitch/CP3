program ising_2d
    implicit none

    integer, parameter :: L = 20
    integer, parameter :: N = L * L
    integer, parameter :: MC_STEPS = 100000

    integer :: i, j, step
    integer :: ip, im, jp, jm
    integer :: dE
    integer :: spin(L, L)

    real(8) :: J, T, beta
    real(8) :: E, M
    real(8) :: r

    call random_seed()

    ! Parameters
    J = 1.0d0
    T = 2.3d0
    beta = 1.0d0 / T

    ! Initialize spins randomly
    do i = 1, L
        do j = 1, L
            call random_number(r)
            if (r < 0.5d0) then
                spin(i,j) = -1
            else
                spin(i,j) = 1
            end if
        end do
    end do

    ! Monte Carlo simulation
    do step = 1, MC_STEPS
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
    E = 0.0d0
    M = 0.0d0

    do i = 1, L
        do j = 1, L
            ip = mod(i, L) + 1
            jp = mod(j, L) + 1

            E = E - J * spin(i,j) * (spin(ip,j) + spin(i,jp))
            M = M + spin(i,j)
        end do
    end do

    print *, "2D Ising Model Simulation"
    print *, "Lattice size =", L, "x", L
    print *, "Temperature =", T
    print *, "Energy per spin =", E / N
    print *, "Magnetization per spin =", M / N

end program ising_2d
