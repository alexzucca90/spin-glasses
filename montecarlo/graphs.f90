module Graphs
    !---------------------------------------------------------------
    !   This module contains some algorithms to generate
    !   the connectivity matrix of some graphs
    !
    !   author: Alex Zucca
    !
    !---------------------------------------------------------------
    use Precision

    implicit none

    contains

    !---------------------------------------------------------------
    !> this subroutine generates the connectivity matrix of
    subroutine GenerateSquare2DLattice(L, Jconn)
        implicit none
            integer, intent(in) :: L
            integer, intent(out) :: Jconn(L**2, L**2)

            integer :: i, j

            !> initialize connection matrix
            Jconn = 0

            do i = 1, L

                if ( i == 1 ) then      !< first line
                    !> upper left corner
                    Jconn(1,L+1) = 1
                    Jconn(1,2) = 1

                    !> upper right corner
                    Jconn(L,2*L) = 1
                    Jconn(L,L-1) = 1

                    !> inner first line
                    do j = 2, L-1
                        Jconn(j, j-1)   = 1
                        Jconn(j, j+1)   = 1
                        Jconn(j, j+L)   = 1
                    end do

                else if ( i == L ) then !< last line
                    !> lower left corner
                    Jconn((L*(L-1)+1),(L*(L-1)+2)) = 1
                    Jconn((L*(L-1)+1),(L*(L-2)+1)) = 1

                    !> lower right corner
                    Jconn(L*L ,L*L-1) = 1
                    Jconn(L*L, L*(L-1)) = 1

                    !> inner last line
                    do j=2, L-1
                        Jconn(L*(L-1)+j, L*(L-1)+j-1 ) = 1
                        Jconn(L*(L-1)+j, L*(L-1)+j+1 ) = 1
                        Jconn(L*(L-1)+j, L*(L-1)+j-L ) = 1
                    end do

                else                    !< bulk
                    !> left side
                    Jconn((i-1)*L+1, (i-1)*L+2) = 1
                    Jconn((i-1)*L+1, (i-2)*L+1) = 1
                    Jconn((i-1)*L+1, i*L+1) = 1

                    !> right side
                    Jconn(i*L, i*L-1) = 1
                    Jconn(i*L, (i-1)*L) = 1
                    Jconn(i*L, (i+1)*L) = 1

                    !> inner line
                    do j = 2, L-1
                        Jconn((i-1)*L+j, (i-1)*L+j-1) = 1
                        Jconn((i-1)*L+j, (i-1)*L+j+1) = 1
                        Jconn((i-1)*L+j, (i-1)*L+j-L) = 1
                        Jconn((i-1)*L+j, (i-1)*L+j+L) = 1
                    end do
                end if

            end do


    end subroutine GenerateSquare2DLattice

    !---------------------------------------------------------------


end module Graphs
