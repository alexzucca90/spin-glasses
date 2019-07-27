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
        integer, intent(in) :: n        !< number of spins in the state
        integer, intent(out) :: s(n)    !< spins array
        integer :: i_s                  !< do loop index
        real(sl) :: rnd                 !< random number

        do i_s = 1, n
            call random_number(rnd)
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
    subroutine MetropolisStep(n, nc, s, J, h, beta, deltaEnergy)

        implicit none

        integer :: n                            !< number of spins
        integer :: nc                           !< number of couplings
        integer, intent(inout) :: s(n)          !< configuration array
        real(dl), intent(in) :: J(nc,3), h(n)   !< coupling and biases
        real(dl), intent(in) :: beta            !< inverse temperature parameter
        real(dl) :: rnd                         !< random number
        real(dl) :: dE                          !< delta Energy
        integer :: i_spin, this_spin            !< spin indexes
        integer :: trial_s(n)                   !< trial configuration

        external deltaEnergy                    !< subroutine to compute the delta energy
        !-------------------------------------------------
        !   the deltaEnergy subroutine is structured in this way
        !
        !   subroutine deltaEnergy(n, nc, s1, s2, J, h, dE)
        !       implicit none
        !       integer, intent(in) :: n
        !       integer, intent(in) :: nc
        !       integer, intent(in) :: s1(n), s2(n)
        !       real(dl), intent(in) :: J(n,n), h(n)
        !       real(dl), intent(out) :: dE
        !
        !        calculate dE
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
            call deltaEnergy(n, nc, trial_s, s, J, h, dE)

            !> accept or reject the update
            call random_number(rnd)
            if ( exp(-dE*beta) .ge. rnd ) then
                s = trial_s
            end if

        end do

    end subroutine MetropolisStep

    !---------------------------------------------------------------
    !> this subroutine swap configurations in Parallel Tempering algorithms
    subroutine SwapConfigurations(n, nc, s1, s2, beta1, beta2, J, h, deltaEnergy)
        implicit none

        integer, intent(in) :: n                !< number of spins
        integer, intent(in) :: nc               !< number of couplings
        integer, intent(inout) :: s1(n), s2(n)  !< spin arrays 1 and 2
        real(dl), intent(in) :: beta1, beta2    !< inverse temperatures 1 and 2
        real(dl), intent(in) :: J(nc,3), h(n)   !< couplings and biases

        external :: deltaEnergy                 !< subroutine that calculates the difference in energy of the two configurations

        real(dl) :: dE
        real(sl) :: rnd
        integer :: stemp(n)

        call deltaEnergy(n, nc, s1, s2, J, h, dE)

        call random_number(rnd)

        if ( exp((beta1 - beta2)*dE) .ge. rnd ) then
            !> swap configurations
            stemp = s1
            s1 = s2
            s2 = stemp
        end if

    end subroutine SwapConfigurations


    !---------------------------------------------------------------
    !> this subroutine implements one instance of simulated annealing
    subroutine SimulatedAnnealing(n, nc, h, J, n_beta, betas, n_equil, s, deltaEnergy)
        implicit none
        integer, intent(in) :: n_beta           !< number of annealing points
        real(dl), intent(in) :: betas(n_beta)   !< annealing schedule
        integer, intent(in) :: n_equil          !< number of equilibration steps
        integer, intent(in) :: n                !< size of the spin-glass system
        integer, intent(in) :: nc               !< number of couplings
        real(dl), intent(in) :: J(nc,3), h(n)   !< coupling and biases of the system
        integer, intent(out) :: s(n)            !< spin configuration

        external :: deltaEnergy

        integer :: i_beta
        integer :: i_equil
        integer :: i_s

        real(dl) :: rnd
        integer :: i
        real(dl) :: beta

        !> randomly initialize the state
        call InitializeSpins(n, s)


        !> follow the annealing schedule
        do i_beta=1, n_beta

            !> select the inverse temperature
            beta = betas( i_beta )

            !> equilibrate the system
            do i_equil=1, n_equil
                !> update the configuration according to the Metropolis-Hastings algorithm
                call MetropolisStep( n, nc, s, J, h, beta, deltaEnergy )
            end do

        end do


    end subroutine SimulatedAnnealing


    !---------------------------------------------------------------
    !> this subroutine calls the simulated annealing instance to fill a sampling array
    ! note: it's not doing the annealing importance sampling
    subroutine SimulatedAnnealingSampling(n, nc, h, J, n_beta, betas, n_equil, n_samples, spins, deltaEnergy)
        implicit none
        integer, intent(in) ::  n_beta
        real(dl), intent(in) :: betas(n_beta)
        integer, intent(in) ::  n_equil
        integer, intent(in) :: n
        integer, intent(in) :: nc
        real(dl), intent(in) :: J(nc,3), h(n)
        integer, intent(in) :: n_samples
        integer, intent(out) :: spins(n_samples, n)

        integer :: spin(n)

        external :: deltaEnergy

        integer :: i_sample

        !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC), PRIVATE(i_sample)
        do i_sample = 1, n_samples
            !> call the simulated annealing instance
            call SimulatedAnnealing(n, nc, h, J, n_beta, betas, n_equil, spins(i_sample,:), deltaEnergy)
        end do
        !$OMP END PARAllEl DO

    end subroutine SimulatedAnnealingSampling


    !---------------------------------------------------------------
    subroutine ParallelTempering(n, nc, h, J, n_replicas, betas, n_sweeps, n_burn_in, n_equil, s, deltaEnergy)
        implicit none

        integer, intent(in) :: n                    !< number of spins
        integer, intent(in) :: nc                   !< number of couplings
        integer, intent(in) :: n_replicas           !< number of replicas
        real(dl), intent(in) :: betas(n_replicas)   !< temperatures of each replica (last one with smaller temperature)
        integer, intent(in) :: n_sweeps             !< number of PT sweeps
        integer, intent(in) :: n_equil              !< number of steps to equilibrate
        real(dl), intent(in) :: J(nc,3)              !< couplings matrix
        real(dl), intent(in) :: h(n)                !< biases matrix
        integer, intent(out) :: s(n_sweeps, n)      !< spins solution
        integer, intent(in) :: n_burn_in            !< number of burn in steps

        external :: deltaEnergy                     !< function that computes the energy

        integer :: spins(n_replicas, n)             !< array of all spins
        integer :: s1(n), s2(n)                     !< spin arrays for single replicas
        integer :: i, k                             !< indexes for do loops
        integer :: i_equil                          !< index for equilibration loop
        real(dl) :: beta1, beta2                    !< inverse temperatures of replicas

        !> begin parallel tempering algorithm
        !> initialize the spins
        do i = 1, n_replicas
            call InitializeSpins(n, s1)
            spins(i,:) = s1
        end do

        do i = 1, n_sweeps+n_burn_in

            !> first equilibrate with Metropolis Hastings for each replica
            !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC), PRIVATE(k)
            do k = 1, n_replicas    !< this can be parallilized
                !s1 = spins(k,:)
                !beta1 = betas(k)
                do i_equil=1, n_equil
                    !> update the configuration according to the Metropolis-Hastings algorithm
                    call MetropolisStep( n, nc, spins(k,:), J, h, betas(k), deltaEnergy )
                end do
                !spins(k,:) = s1
            end do
            !$OMP END PARAllEl DO

            !> swapping configurations
            do k = n_replicas, 2, -1

                !> set up the variables
                beta1 = betas(k)
                beta2 = betas(k-1)
                s1 = spins(k,:)
                s2 = spins(k-1,:)

                call SwapConfigurations(n, nc, s1, s2, beta1, beta2, J, h, deltaEnergy)

                spins(k,:) = s1
                spins(k-1,:) = s2

            end do

            !> attach the solutions if past burn in
            if ( i > n_burn_in ) then
                s(i-n_burn_in,:) = spins(n_replicas,:)
            end if

        end do

        !> solutions (so far, only one solution)
        !s = spins(n_replicas,:)

    end subroutine ParallelTempering


end module MonteCarlo
