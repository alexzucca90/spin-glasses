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
    function Energy(n, s, h, J)

        implicit none

        integer, intent(in) :: n        !< number of spins
        real(dl), intent(in) :: J(n,n)   !< couplings
        real(dl), intent(in) :: h(n)    !< biases
        integer, intent(in) :: s(n)     !< spin configuration
        real(dl) :: Energy              !< energy

        integer :: i, k

        Energy = 0.d0
        do i = 1, n
            Energy = Energy - h(i)*s(i)
            do k=i, n
                Energy = Energy -  J(i,k)*s(i)*s(k)
            end do
        end do
    end function Energy

    !---------------------------------------------------------------
    !> this subroutine calculates the energy difference between two
    !   spin configurations
    subroutine delta_Energy(n, s1, s2, J, h, dE)

        implicit none

        integer, intent(in) :: n
        integer, intent(in) ::  s1(n), s2(n)
        real(dl), intent(in) :: J(n,n)
        real(dl), intent(in) :: h(n)
        real(dl), intent(out) :: dE

        dE = Energy(n, s1, h, J) - Energy(n, s2, h, J)

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
