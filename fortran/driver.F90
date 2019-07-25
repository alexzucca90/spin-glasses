program driver
!---------------------------------------------------------------
!   Simple program that solves a Hen Frustrated Cluster Loop
!   spin-glass Hamiltonian with simulated annealing.
!
!
!---------------------------------------------------------------

    use Graphs
    use Problems
    use Precision
    use MonteCarlo
    use Model

    implicit none

    integer :: L                                        !< size of the problem
    real(dl), dimension(:,:), allocatable :: Jproblem   !< coupling  matrix of the lattice
    integer, dimension(:,:), allocatable :: Jconn       !< graph connectivity matrix
    real(dl), dimension(:), allocatable :: h            !< biases of the model
    integer :: n                                        !< number of spins
    real(dl) :: alpha, R                                !< parameters of the H-FCL
    integer, dimension(:), allocatable :: solution      !< spin configuration of the solution provided
    real(dl), dimension(:), allocatable :: betas        !< array containing the inverse temperatures of the annealing schedule
    integer :: n_beta                                   !< number of temperatures in the schedule
    real(dL) :: exponent, T_max                         !< exponent of the temperature
    integer :: n_equil                                  !< number of steps to equilibrate at each temperature
    integer :: n_samples                                !< number of samples of low-energy states
    integer, dimension(:,:), allocatable :: solutions   !< solution samples
    integer :: n_burn_in                                !< burn in for PT

    integer :: i

    real    :: start, finish        !< for benchmarking the code
    integer :: ThreadNum = 0        !< OMP threads parameter: leave to 0 for full parallel

    !> set the parameters of the H-FCL
    alpha   = 0.8
    R       = 20

    !> set the annealing schedule parameters
    n_beta = 25
    allocate(betas(n_beta))
    T_max       = 10.d0
    exponent    = 0.75d0
    n_equil     = 100

    !> generate the connectivity graph of a Chimera lattice
    L = 3
    n = 8*(L**2)
    allocate(Jconn(n,n))

    call GenerateSquareChimeraGraph(L, Jconn)

    !> writing the connectivity matrix on file
    open(unit=4, file = './problems/chimera_lattice.dat', status = 'unknown', recl=9999)
    do i = 1,  n
        write(4,*) Jconn(i,:)
    end do
    close(4)
    !deallocate(Jconn)

    !> set the number of spins
    !L = 8
    !n = L**2

    !allocate(Jconn(n,n))
    allocate(Jproblem(n, n))
    allocate(h(n))
    allocate(solution(n))

    h = 0
    solution  = 0
    n_samples = 20

    allocate(solutions(n_samples, n))

    !write(*,*) 'Generating 2D lattice.'
    !call GenerateSquare2DLattice(L, Jconn)
    !open(unit=4, file = './problems/square_lattice.dat', status = 'unknown', recl=9999)
    !do i = 1,  n
    !    write(4,*) Jconn(i,:)
    !end do
    !close(4)

    !> generate the annealing schedule
    write(*,*) 'Making annealing schedule'
    call MakeExponentialTemperatureSheduleMaxMin( n_beta, betas, exponent, T_max )

    !> Now generate the problem
    write(*,*) 'Making H-FCL. Hamiltonian:'
    open(unit=1, file = './problems/HFCL_chimera.dat', status='unknown', recl = 9999)
    call HenFrustratedClusterLoop( n, alpha, R, Jconn, Jproblem )
    do i = 1,  n
        write(1,*) Jproblem(i,:)
    end do
    close(1)
    write(*,*)

    !> Now find the candidate solution
    write(*,*) 'Starting simulation'
    !call SimulatedAnnealing( n_beta, betas, n_equil, Jproblem, h, n, solution, delta_Energy )

    !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)

    call cpu_time(start)
    call SimulatedAnnealingSampling( n_beta, betas, n_equil, Jproblem, h, n, n_samples, solutions, delta_Energy)
    call cpu_time(finish)
    print '("Time SA = ",f6.3," seconds.")',finish-start
    !write(*,*) 'writing solutions'
    open(unit=2, file='./output/solutions_SA.dat', status='unknown', recl = 9999)
    do i = 1, n_samples
        !write(*,*) solutions(i,:)
        write(2,*) solutions(i,:)
    end do
    close(2)

    open(unit=3, file = './output/energies_SA.dat', status='unknown')
    do i = 1, n_samples
        !solution = solutions(i,:)
        write(3,*) Energy(n, solutions(i,:), h, Jproblem)
    end do
    close(3)

    !> here add parallel tempering

    solutions = 0
    n_burn_in = 200

    call cpu_time(start)
    call ParallelTempering(n, n_beta, betas, n_samples, n_burn_in, n_equil, Jproblem, h, solutions, delta_Energy)
    call cpu_time(finish)
    print '("Time PT = ",f6.3," seconds.")',finish-start

    open(unit=2, file='./output/solutions_PT.dat', status='unknown', recl = 9999)
    do i = 1, n_samples
    !write(*,*) solutions(i,:)
    write(2,*) solutions(i,:)
    end do
    close(2)

    open(unit=3, file = './output/energies_PT.dat', status='unknown')
    do i = 1, n_samples
    !solution = solutions(i,:)
    write(3,*) Energy(n, solutions(i,:), h, Jproblem)
    end do
    close(3)

end program driver
