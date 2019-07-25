module metric

    use Precision

    implicit none

    integer, parameter :: dl = KIND(1.d0)

contains

    function Energy(n,h,J,s)
        implicit none
        integer :: n
        real(dl) :: Energy
        real(dl) :: h(n), J(n,n)
        integer :: s(n)

        integer :: i,k

        Energy = 0.d0
        do i = 1, n
            Energy = Energy + h(i)*s(i)
        end do

        do i = 1,n-1
            do k = i+1,n
                Energy = Energy + s(i)*s(k)*J(i,k)
            end do
        end do

    end function Energy

    function HammingDistance(n, s1, s2)
        implicit none
        integer :: n
        integer :: s1(n), s2(n)
        integer :: HammingDistance

        integer :: i,j

        HammingDistance = 0
        do i = 1, n
            if (s1(i) /= s2(i)) HammingDistance = HammingDistance + 1
        end do

    end function


end module metric

!----------------------------------------------------------------------
!> This program generate a random problem and computes all the energies
program main

    use metric
    use Graphs
    use problems

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

    integer :: i

    !> first, let's generate the graph connectivity, for example a square lattice with 4 spins per size
    L = 4
    n = L**2

    allocate(Jconn(n,n))

    !> generate the square lattice
    call GenerateSquare2DLattice(L, Jconn)

    !> Now generate a random H-FCL problem
    allocate(Jproblem(n, n))
    allocate(h(n))
    allocate(solution(n))

    h = 0
    solution  = 0

    !> set the H-FCL parameters
    alpha = 0.3
    R = 10
    call HenFrustratedClusterLoop( n, alpha, R, Jconn, Jproblem )

    !> checking how many spins are active





    real(dl), dimension( n_spins, n_spins ) :: J
    real(dl), dimension( n_spins ) :: h
    integer, dimension( n_spins ) :: s

    real(dl) ::  alpha = 0.1
    real(dl) :: R = 0.9
    real(dl) :: remainder

    real(dl), dimension( 2**n_spins ) :: energies
    integer :: i, k, v, w, z, m

    integer, dimension(2**n_spins,n_spins) :: all_states
    integer, dimension(:,:), allocatable :: fit_samples
    integer, dimension(:,:), allocatable :: adjacency_matrix

    real(dl) :: min_en, max_en, energy_threshold
    real(dl) :: dist


    J = 0

    h = (/1, -1, 0, 1, -1, -1, -1, -1, 1, 0, -1, 1, 1, 1, 1, 1/)
    J(1,2) = -1
    J(2,1)  = -1
    J(2,4) = -1
    J(2,5) = -1
    J(2,6) = -1
    J(3,4) = -1
    J(3,7) = -1
    J(3,8) = -1
    J(5,10) = -1
    J(5,13) = -1
    J(6,10) = -1
    J(6,14) = -1
    J(7,15) = -1
    J(8,10) = -1
    J(8,16) = -1
    J(9,10) = -1
    J(10,11) = -1
    J(11,12) = -1

    write(*,*) ""

    ! open a file to store the exact search
    open(unit = 1, file = 'energies.dat', status = 'unknown', recl = 9999)

    ! First, check which states have the energy below the threshold
    do i = 1, 2**n_spins

        remainder = i
        do k = n_spins, 1, -1

            s(k) = int(remainder/2**(k-1))
            remainder = mod(int(remainder), 2**(k-1))
            !write(*,*) s(k), remainder
            if ( s(k) < 1 ) then
                s(k) = -1
            else if ( s(k) > 1 ) then
                s(k) = -1
            end if

            all_states(i,:) = s

        end do


        energies(i) = Energy(n_spins, h, J, s)
        write(1,*) i, energies(i), s

    end do
    close(1)

    ! Now find the maximum and minimum energies
    min_en = MINVAL(energies)
    max_en = MAXVAL(energies)

    write(*,*) min_en, max_en

    energy_threshold = min_en + alpha*(max_en - min_en)


    open(unit=1, file = "fit_samples.dat", status="unknown", recl=99999)
    k = 0
    do i = 1, 2**n_spins
        if ( energies(i) .le. energy_threshold ) then
            !write(*,*) i
            write(1,*) all_states(i,:)
            k = k+1
        end if
    end do
    close(1)

    write(*,*) "Energy threshold:", energy_threshold, "number of fit samples:", k

    allocate(fit_samples(k,n_spins))
    allocate(adjacency_matrix(k,k))

    adjacency_matrix = 0

    m=1
    do i = 1, 2**n_spins
        if ( energies(i) .le. energy_threshold ) then
            fit_samples(m,:) = all_states(i,:)
            m = m+1
        end if
    end do


    do i=1,k-1
        do m=i,k
            dist = HammingDistance(n_spins, fit_samples(i,:), fit_samples(m,:))
            if(dist > n_spins*R) then
                adjacency_matrix(i,m) = 1
            end if
        end do
    end do

    ! write the adjacency_matrix on file
    open(unit=1, file = "adjacency_matrix.dat", status = "unknown", recl = 99999)
    open(unit=2, file = "edges_list.dat", status="unknown")
    do i = 1, k
        write(1,*) adjacency_matrix(i,:)
        do m = i,k
            if ( adjacency_matrix(i,m) > 0 ) then
                write(2,*) i, m
            end if
        end do
    end do
    close(1)
    close(2)

    write(*,*) "Number of connections:", SUM(adjacency_matrix)

end program main
