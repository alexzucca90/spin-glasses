module MonteCarlo
!---------------------------------------------------------------
!   This module contains some Monte Carlo subroutines to
!   solve and sample spin-glass systems
!
!   author: Alex Zucca
!
!---------------------------------------------------------------

    use Precision

    implicit none

    contains

    !---------------------------------------------------------------
    !> this subroutine initializes the spins randomly up or down
    subroutine InitializeSpins(n, s)
        implicit none
        integer, intent(in) :: n    !< number of spins in the state
        integer, intent (out) :: s  !< spins array
        integer :: i_s              !< do loop index
        real(sl) :: rnd             !< random number

        do i_s = 1, n
            write(*,*) rnd
            if ( rnd > 0.5 ) then
                s( i_s ) = 1
            else
                s( i_s ) = -1
            end if
        end do

    end subroutine InitializeSpins


    !---------------------------------------------------------------
    !> this subroutine performs a Metropolis-Hastings update of the
    !  spin glass system
    subroutine MetropolisStep(n, s, J, h, beta, deltaEnergy)

        implicit none

        integer :: n                            !< number of spins
        integer, intent(inout) :: s(n)          !< configuration array
        real(dl), intent(in) :: J(n,n), h(n)    !< coupling and biases
        real(dl), intent(in) :: beta            !< inverse temperature parameter
        real(dl) :: rnd                         !< random number
        real(dl) :: dE                          !< delta Energy
        integer :: i_spin, this_spin            !< spin indexes
        integer :: trial_s(n)                   !< trial configuration

        external deltaEnergy                    !< subroutine to compute the delta energy
        !-------------------------------------------------
        !   subroutine deltaEnergy(trial_s, s, J, h, dE)
        !       implicit none
        !
        !
        !   end subroutine deltaEnergy
        !-------------------------------------------------

        !> perform n updates for each Metropolis-Hastings step
        do i_spin = 1, n

            trial_s = s

            ! randomly select the spin to flip
            call random_number(rnd)
            this_spin = int(n*rnd) + 1

            !> flip the spin
            trial_s(this_spin) = - s(this_spin)

            !> calculate the delta of the energy configuration
            call deltaEnergy(n, trial_s, s, J, h, dE)

            !> accept or reject the update
            call random_number(rnd)
            if ( exp(-dE*beta) .ge. rnd ) then
                s = trial_s
            end if

        end do

    end subroutine MetropolisStep


    !---------------------------------------------------------------
    !> this subroutine implements one instance of simulated annealing
    subroutine SimulatedAnnealing(n_beta, betas, n_equil, J, h, n, s, deltaEnergy)
        implicit none
        integer, intent(in) :: n_beta           !< number of annealing points
        real(dl), intent(in) :: betas(n_beta)    !< annealing schedule
        integer, intent(in) :: n_equil             !< number of equilibration steps
        integer, intent(in) :: n                !< size of the spin-glass system
        real(dl), intent(in) :: J(n,n), h(n)    !< coupling and biases of the system
        integer, intent(out) :: s(n)            !< spin configuration

        external :: deltaEnergy

        integer :: i_beta
        integer :: i_equil
        integer :: i_s

        real(dl) :: rnd
        integer :: i
        real(dl) :: beta

        !write(*,*) 'Setting initial random state:'
        call InitializeSpins(n, s)


        !> follow the annealing schedule
        !write(*,*) 'Starting annealing:'
        do i_beta=1, n_beta

            !> select the inverse temperature
            beta = betas( i_beta )

            !write(*,*) '     temperature:', 1.d0/beta

            !> equilibrate the system
            do i_equil=1, n_equil
                !> update the configuration according to the Metropolis-Hastings algorithm
                call MetropolisStep( n, s, J, h, beta, deltaEnergy )
            end do

        end do

    end subroutine SimulatedAnnealing


    !---------------------------------------------------------------
    !> this subroutine calls the simulated annealing instance to fill a sampling array
    ! note: it's not doing the annealing importance sampling
    subroutine SimulatedAnnealingSampling(n_beta, betas, n_equil, J, h, n, n_samples, spins, deltaEnergy)
        implicit none
        integer, intent(in) ::  n_beta
        real(dl), intent(in) :: betas(n_beta)
        integer, intent(in) ::  n_equil
        integer, intent(in) :: n
        real(dl), intent(in) :: J(n,n), h(n)
        integer, intent(in) :: n_samples
        integer, intent(out) :: spins(n_samples, n)

        integer :: spin(n)

        external :: deltaEnergy

        integer :: i_sample

        do i_sample = 1, n_samples
            !> call the simulated annealing instance
            call SimulatedAnnealing(n_beta, betas, n_equil, J, h, n, spin, deltaEnergy)
            spins(i_sample,:) = spin
        end do

    end subroutine SimulatedAnnealingSampling


    !---------------------------------------------------------------
    subroutine ParallelTempering(n_replicas, betas, n_sweeps, n_equil, J, h, n, s, deltaEnergy)
        implicit none

        integer, intent(in) :: n_replicas           !< number of replicas
        real(dl), intent(in) :: betas(n_replicas)   !< temperatures of each replica (last one with smaller temperature)
        integer, intent(in) :: n_sweeps             !< number of PT sweeps
        integer, intent(in) :: n_equil              !< number of steps to equilibrate
        real(dl), intent(in) :: J(n,n)              !< couplings matrix
        real(dl), intent(in) :: h(n)                !< biases matrix
        real(dl), intent(out) :: s(n)               !< spins solution

        external :: deltaEnergy                     !< function that computes the energy

        integer :: spins(n_replicas, n)             !<
        integer :: spin_replica(n)
        integer :: i, j
        integer :: i_equil
        real(dl) :: beta

        !> begin parallel tempering algorithm
        !> initialize the spins
        do i = 1, n_replicas
            call InitializeSpins(n, spin_replica)
            spins(i,:) = spin_replica
        end do

        do i = 1, n_sweeps
            !> first equilibrate with Metropolis Hastings for each replica
            do j = 1, n_replicas    !< this can be parallilized
                spin_replica = spins(j,:)
                do i_equil=1, n_equil
                    !> update the configuration according to the Metropolis-Hastings algorithm
                    beta = betas
                    call MetropolisStep( n, spin_replica, J, h, beta, deltaEnergy )
                end do
                spins(j,:) = spin_replica
            end do

            !> swapping configurations
            do j = 2, n_replicas

            end do

        end do



    end subroutine ParallelTempering


end module MonteCarlo
