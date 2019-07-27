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
            integer, intent(out) :: Jconn(L*L, L*L)

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
    !> This subroutine generates the connectivity matrix of a square Chimera graph
    subroutine GenerateSquareChimeraGraph(L, Jconn)
        implicit none
        integer, intent(in) :: L                        !< number of unit cells on the side
        integer, intent(out) :: Jconn(8*L**2, 8*L**2)   !< connectivity matrix

        integer :: i, j, k

        !> initialize connection matrix
        Jconn = 0

        !> start with the intercell connections
        do i = 1, L
            !> first line
            if ( i == 1 ) then

                do k = 1, 4
                    Jconn(k,8*L+k) = 1                  !< upper left, vertical
                    Jconn(4+k, 12+k) = 1                !< upper left, horizontal
                    Jconn(8*L-(8-k), 16*L-(8-k)) = 1    !< upper right, vertical
                    Jconn(8*L-(k-1), 8*(L-1)-(k-1)) = 1 !< upper right, horizontal
                end do


                do j = 2,L-1
                    do k = 1, 4
                        Jconn(8*(j-1)+k, 8*(L+j-1)+k) = 1   !< vertical, down
                        Jconn(8*(j-1)+4+k, 8*(j-2)+4+k) = 1 !< horizontal, left
                        Jconn(8*(j-1)+4+k, 8*j+4+k) = 1     !< horizontal, right
                    end do
                end do

            !> last line
            else if ( i == L ) then

                do k = 1, 4
                    Jconn(8*L*(L-1)+k, 8*L*(L-2)+k ) = 1        !< lower left, vertical
                    Jconn(8*L*(L-1)+4+k, 8*l*(L-1)+12+k) = 1    !< lower left, horizontal
                    Jconn(8*L**2-(8-k), 8*L*(L-1)-(8-k)) = 1    !< lower right, veritical
                    Jconn(8*L**2-(k-1), 8*L**2-8-(k-1)) = 1     !< lower right, horizontal
                end do

                do j = 2, L-1
                    do k = 1, 4
                        Jconn(8*L*(L-1)+8*(j-1)+k, 8*L*(L-2)+8*(j-1)+k) = 1     !< vertical up
                        Jconn(8*L*(L-1)+8*(j-1)+4+k, 8*L*(L-1)+8*(j-2)+4+k) = 1 !< horizontal, left
                        Jconn(8*L*(L-1)+8*(j-1)+4+k, 8*L*(L-1)+8*(j)+4+k) = 1   !< horizontal, right
                    end do
                end do

            !> bulk
            else

                do k = 1, 4
                    Jconn(8*L*(i-1)+k, 8*L*(i-2)+k) = 1     !< left side, vertical up
                    Jconn(8*L*(i-1)+k, 8*L*(i)+k) = 1       !< left side, vertical down
                    Jconn(8*L*(i-1)+4+k,8*L*(i-1)+12+k) = 1 !< left side, horizontal right
                    Jconn(8*L*i-(8-k), 8*L*(i-1)-(8-k)) = 1 !< right side, vertical up
                    Jconn(8*L*i-(8-k), 8*L*(i+1)-(8-k)) = 1 !< right side, vertical down
                    Jconn(8*L*i-(k-1), 8*L*i-8-(k-1)) = 1   !< right side, horizontal left
                end do

                do j = 2,L-1
                    do k = 1, 4
                        Jconn(8*L*(i-1)+8*(j-1)+k, 8*L*(i-2)+8*(j-1)+k) = 1     !< vertical, up
                        Jconn(8*L*(i-1)+8*(j-1)+k, 8*L*(i)+8*(j-1)+k) = 1       !< vertical, down
                        Jconn(8*L*(i-1)+8*(j-1)+4+k, 8*L*(i-1)+8*(j-2)+4+k) = 1 !< horizontal, left
                        Jconn(8*L*(i-1)+8*(j-1)+4+k, 8*L*(i-1)+8*(j)+4+k) = 1   !< horizontal, right
                    end do
                end do

            end if

        end do

        !> Now add the intra-cell connections
        do i = 1, L**2
            do j=1,4
                do k=1,4
                    Jconn(8*(i-1)+j, 8*(i-1)+4+k) = 1   !< horizontal-vertical
                    Jconn(8*(i-1)+4+k, 8*(i-1)+j) = 1   !< vertical-horizontal
                end do
            end do
        end do

    end subroutine GenerateSquareChimeraGraph

end module Graphs
