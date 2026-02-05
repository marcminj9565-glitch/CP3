program ising_1d_energy_heat
    implicit none

    integer :: N
    real(8) :: J, T, beta, kB
    real(8) :: U_pbc, U_obc
    real(8) :: C_pbc, C_obc
    real(8) :: sech2

    ! Parameters
    J  = 1.0d0
    kB = 1.0d0
    T  = 2.0d0
    beta = 1.0d0 / (kB * T)

    ! Large system size
    N = 1000000

    ! sech^2(beta J)
    sech2 = 1.0d0 / (cosh(beta * J)**2)

    ! Internal energy per spin
    U_pbc = -J * tanh(beta * J)
    U_obc = -J * (1.0d0 - 1.0d0 / N) * tanh(beta * J)

    ! Heat capacity per spin
    C_pbc = kB * (beta * J)**2 * sech2
    C_obc = kB * (1.0d0 - 1.0d0 / N) * (beta * J)**2 * sech2

    print *, "1D Ising Model (N =", N, ")"
    print *, "Temperature T =", T
    print *, "-------------------------------------"
    print *, "Internal Energy per spin:"
    print *, "PBC :", U_pbc
    print *, "OBC :", U_obc
    print *, "-------------------------------------"
    print *, "Heat Capacity per spin:"
    print *, "PBC :", C_pbc
    print *, "OBC :", C_obc
    print *, "-------------------------------------"
    print *, "Difference (PBC - OBC):"
    print *, "ΔU =", U_pbc - U_obc
    print *, "ΔC =", C_pbc - C_obc

end program ising_1d_energy_heat
