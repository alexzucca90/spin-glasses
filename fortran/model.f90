module Model
    !---------------------------------------------------------------
    !   This module contains the model definition for the spin-glass
    !   system, and some annealing schedules generators.
    !
    !   author: Alex Zucca
    !
    !---------------------------------------------------------------
    use Precision

    implicit none

    contains

    !---------------------------------------------------------------
    !> this function calculates the energy for a spin configuration
    !   and a Hamiltonian
    function Energy(n, nc, s, h, J)

        implicit none

        integer, intent(in) :: n        !< number of spins
        integer, intent(in) :: nc       !< number of couplings
        real(dl), intent(in) :: J(nc,3) !< couplings
        real(dl), intent(in) :: h(n)    !< biases
        integer, intent(in) ::  s(n)    !< spin configuration
        real(dl) :: Energy              !< energy

        integer :: i, k, ic !< some indeces
        real(dl) :: Jik     !< coupling value

        Energy = 0.d0
        Energy = Energy + SUM(h*s)
        do ic=1, nc
            i = int(J(ic,1))
            k = int(J(ic,2))
            Jik = J(ic,3)
            Energy = Energy+ Jik*s(i)*s(k)
        end do

    end function Energy



    !---------------------------------------------------------------
    !> this subroutine calculates the energy difference between two
    !   spin configurations
    subroutine delta_Energy(n, nc, s1, s2, J, h, dE)

        implicit none

        integer, intent(in) :: n                !< number of spins
        integer, intent(in) :: nc               !< number of couplings
        integer, intent(in) ::  s1(n), s2(n)    !< spins cofigurations
        real(dl), intent(in) :: J(nc,3)          !< couplings
        real(dl), intent(in) :: h(n)            !< biases
        real(dl), intent(out) :: dE             !< delta_energy

        dE = Energy(n, nc, s1, h, J) - Energy(n, nc, s2, h, J)

    end subroutine

    !---------------------------------------------------------------
    !> this subroutine genereates a linear annealing schedule
    subroutine MakeLinearTemperatureSchedule(n_beta, betas, T_min, T_max)
        implicit none

        integer, intent(in) :: n_beta           !< number of temperatures in the annealing schedule
        real(dl), intent(out) :: betas(n_beta)  !< arrays containing the annealing schedule
        real(dl), intent(in) :: T_min, T_max    !< min and max temperatures (k_b = 1) of the annealing

        integer :: i

        do i = 1, n_beta
            betas(i) = 1.d0/ T_max+ (i-1)*(1.d0/T_min - 1.d0/T_max)/(n_beta-1)
        end do

    end subroutine MakeLinearTemperatureSchedule

    !---------------------------------------------------------------
    !> this subroutines generates an exponential annealing schedule -starting from T max
    subroutine MakeExponentialTemperatureSheduleMaxMin(n_beta, betas, exponent, T_max)
        implicit none
        integer, intent(in) :: n_beta
        real(dl), intent(out) :: betas(n_beta)
        real(dl), intent(in) :: exponent
        real(dl), intent(in) :: T_max

        real(dl) :: this_T
        integer :: i

        if (exponent .ge. 1.d0) then
            write(*,*) 'WARNING: temperature increasing in the schedule'
            stop
        end if

        this_T = T_max
        betas(1) = 1.d0 / this_T

        do i = 2, n_beta
            this_T = this_T**exponent
            betas(i) = 1.d0 / this_T
        end do

    end subroutine MakeExponentialTemperatureSheduleMaxMin

    !---------------------------------------------------------------
    !> this subroutines generates an exponential annealing schedule starting from T_min
    subroutine MakeExponentialTemperatureSheduleMinMax(n_beta, betas, exponent, T_min)
        implicit none
        integer, intent(in) :: n_beta
        real(dl), intent(out) :: betas(n_beta)
        real(dl), intent(in) :: exponent
        real(dl), intent(in) :: T_min

        real(dl) :: this_T
        integer :: i

        if (exponent .ge. 1.d0) then
            write(*,*) 'WARNING: temperature increasing in the schedule'
            stop
        end if

        this_T = T_min
        betas(n_beta) = 1.d0 / this_T

        do i = n_beta - 1, 1, -1
            this_T = this_T**(1.d0/exponent)
            betas(i) = 1.d0 / this_T
        end do

    end subroutine MakeExponentialTemperatureSheduleMinMax


end module Model
