module Problems
    !---------------------------------------------------------------
    !   This module contains some algorithms to prepare some spin-glass
    !   problems
    !
    !   author: Alex Zucca
    !
    !---------------------------------------------------------------
    use Precision

    implicit none

    contains

    !---------------------------------------------------------------
    !> This subroutine generates a la Hen Frustrated Cluster Loop (H-FCL)
    !   problem Hamiltonian used to benchmark classical
    !   and quantum adiabatic algorithms
    subroutine HenFrustratedClusterLoop(n, alpha, R, Jconn, Jproblem)
        implicit none
        integer, intent(in) :: n                !< number of nodes (spins) of the system
        real(dl), intent(in) :: alpha           !< density of constraints for the problem
        real(dl), intent(in) :: R               !< ruggedness of the problem
        integer, intent(in) :: Jconn(n,n)       !< connectivity matrix of the graph (for D Wave it has to be the Chimera or Pegasus)
        real(dl), intent(out) :: Jproblem(n,n)   !< problem Hamiltonian couplings

        !> internal variables
        integer :: M                !< M = alpha* n
        integer :: i_loop           !< index for loop
        integer :: this_spin        !< current spin index
        integer :: previous_spin    !< previous spin index
        integer :: nn_num           !< number of nearest neighbours of the spin
        integer :: i_walker         !< steps done by the random walker
        integer :: loop_path(n)     !< array containing the spins visited by the random walker
        integer :: loop_counter(n)  !< counting how many times the walker passed through the spin
        real(dl) :: rnd             !< random number instance
        integer :: i_nn             !< nearest neighbour index
        integer :: i,j              !< extra index for do loops
        integer :: index            !< index of the loop to add frustration
        integer, dimension(:), allocatable :: nearest_neighbours    !< nearest neighbours array
        integer :: start_index, end_index
        logical :: open_loop

100     M = int(alpha*n)+1

        i_loop = 0
        do while ( i_loop <  M )

            write(*,*) 'creating loop:', i_loop

            !> reset the arrays
            loop_counter    = 0
            loop_path       = 0

            !> select a random starting point
            call random_number(rnd)
            this_spin = int( rnd * n ) + 1
            !> set the previous spin to 0
            previous_spin = 0

            !> adding the spin counters
            loop_path(1) = this_spin
            loop_counter(this_spin) = 1

            !> random walker step
            i_walker = 1

            open_loop = .true.

            do while ( i_walker < n .and. open_loop )

                !> start by checking how many nearest neighbours the spin has
                nn_num = sum(Jconn(this_spin,:))
                if ( previous_spin /= 0 ) then
                    nn_num = nn_num - 1
                end if

                if ( nn_num == 0 ) then
                    write(*,*) 'WARNING: the loop is at a dead end: restarting.'
                    goto 100
                end if

                !> allocate and initialize the nearest neighbours array
                allocate(nearest_neighbours(nn_num))
                nearest_neighbours = 0

                !> fill the nearest neighbours array
                i_nn = 1
                do i = 1, n
                    if ( Jconn(this_spin,i) > 0) then
                        if ( i /= previous_spin ) then
                            nearest_neighbours(i_nn) = i
                            i_nn = i_nn+1
                        end if
                    end if
                end do

                !> select the next spin of the path
                call random_number(rnd)
                previous_spin = this_spin
                this_spin = nearest_neighbours(int(rnd*nn_num) + 1)

                !> de-allocate the nearest neighbours array
                deallocate(nearest_neighbours)

                i_walker = i_walker+1

                !> add the new spin to the path and to the counter
                loop_path(i_walker) = this_spin
                loop_counter(this_spin) = loop_counter(this_spin) + 1


                !> check whether the loop is closed
                if ( loop_counter( this_spin ) > 1 ) then !< the walker returned to some previous point.
                    open_loop = .false.
                    i_loop = i_loop + 1
                end if

            end do

            !> loop is closed, so find the extremes of the loop
            do i = 1, n
                if ( loop_path(i) == this_spin ) then
                    start_index = i
                    do j=i+1, n
                        if ( loop_path(j) == this_spin )  then
                            end_index = j
                            exit
                        end if
                    end do
                    exit
                end if
            end do

            !> now set the couplings of the loops to one
            do i = start_index+1, end_index
                if ( loop_path(i-1) < loop_path(i) ) then
                    Jproblem(loop_path(i-1), loop_path(i)) = Jproblem(loop_path(i-1), loop_path(i)) - 1
                else
                    Jproblem(loop_path(i), loop_path(i-1)) = Jproblem(loop_path(i), loop_path(i-1)) - 1
                end if
            end do

            !> insert the frustration on the loop in a random coupling
            call random_number(rnd)
            index = start_index + int(rnd*(end_index - start_index)) + 1
            if (loop_path(index-1) < loop_path(index)) then
                Jproblem(loop_path(index-1), loop_path(index)) = Jproblem(loop_path(index-1), loop_path(index)) + 2
            else
                Jproblem(loop_path(index), loop_path(index-1)) = Jproblem(loop_path(index-1), loop_path(index)) + 2
            end if

        end do

        !> check whether the loop is acceptable
        do i = 1, n
            do j = i, n
                if ( abs(Jproblem(i,j)) > R ) then
                    write(*,*) 'Loop exceeds ruggedness. Discarding'
                    goto 100
                end if
            end do
        end do

    end subroutine HenFrustratedClusterLoop

    !---------------------------------------------------------------
    !> This subroutine generates a la King Frustrated Cluster Loop (K-FCL)
    !   problem Hamiltonian used to benchmark classical
    !   and quantum adiabatic algorithms. K-FCL are used only on the Chimera graph
    subroutine KingFrustratedClusterLoop
        implicit none
    end subroutine KingFrustratedClusterLoop


    !---------------------------------------------------------------
    !> This subroutine generates a Deceptive Cluster Loop (H-FCL)
    !   problem Hamiltonian used to benchmark classical
    !   and quantum adiabatic algorithms
    subroutine DeceptiveClusterLoop
        implicit none
    end subroutine DeceptiveClusterLoop


end module Problems