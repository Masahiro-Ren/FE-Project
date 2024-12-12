!OCL SERIAL
subroutine apply_filter_xyz_direction_orig(filterMat, filterMat_tr, filterMat_vtr, q_in, q_tmp )
    implicit none

    real(8), intent(in) :: filterMat(12, 12)
    real(8), intent(in) :: filterMat_tr(12, 12)
    real(8), intent(in) :: filterMat_vtr(12, 12)
    real(8), intent(inout) :: q_in(12, 12, 12)
    real(8), intent(inout) :: q_tmp(12, 12, 12)
    
    integer :: i, j, k

    !-- x direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

    q_tmp(i,j,k) = filterMat(i,1)  * q_in(1,j,k) + &
                    filterMat(i,2)  * q_in(2,j,k) + &
                    filterMat(i,3)  * q_in(3,j,k) + & 
                    filterMat(i,4)  * q_in(4,j,k) + &
                    filterMat(i,5)  * q_in(5,j,k) + & 
                    filterMat(i,6)  * q_in(6,j,k) + & 
                    filterMat(i,7)  * q_in(7,j,k) + & 
                    filterMat(i,8)  * q_in(8,j,k) + &
                    filterMat(i,9)  * q_in(9,j,k) + & 
                    filterMat(i,10) * q_in(10,j,k) + & 
                    filterMat(i,11) * q_in(11,j,k) + & 
                    filterMat(i,12) * q_in(12,j,k)  

    end do
    end do
    end do

    !-- y direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

    q_in(i,j,k) = q_tmp(i,1,k)  * filterMat_tr(1,j) + &
                    q_tmp(i,2,k)  * filterMat_tr(2,j) + &
                    q_tmp(i,3,k)  * filterMat_tr(3,j) + &
                    q_tmp(i,4,k)  * filterMat_tr(4,j) + &
                    q_tmp(i,5,k)  * filterMat_tr(5,j) + &
                    q_tmp(i,6,k)  * filterMat_tr(6,j) + &
                    q_tmp(i,7,k)  * filterMat_tr(7,j) + &
                    q_tmp(i,8,k)  * filterMat_tr(8,j) + &
                    q_tmp(i,9,k)  * filterMat_tr(9,j) + &
                    q_tmp(i,10,k) * filterMat_tr(10,j) + &
                    q_tmp(i,11,k) * filterMat_tr(11,j) + &
                    q_tmp(i,12,k) * filterMat_tr(12,j) 

    end do
    end do
    end do

    !-- z direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

    q_tmp(i,j,k) = q_in(i,j,1)  * filterMat_vtr(1,k) + &
                    q_in(i,j,2)  * filterMat_vtr(2,k) + &
                    q_in(i,j,3)  * filterMat_vtr(3,k) + &
                    q_in(i,j,4)  * filterMat_vtr(4,k) + &
                    q_in(i,j,5)  * filterMat_vtr(5,k) + &
                    q_in(i,j,6)  * filterMat_vtr(6,k) + &
                    q_in(i,j,7)  * filterMat_vtr(7,k) + &
                    q_in(i,j,8)  * filterMat_vtr(8,k) + &
                    q_in(i,j,9)  * filterMat_vtr(9,k) + &
                    q_in(i,j,10) * filterMat_vtr(10,k) + &
                    q_in(i,j,11) * filterMat_vtr(11,k) + &
                    q_in(i,j,12) * filterMat_vtr(12,k) 

    end do
    end do
    end do
end subroutine apply_filter_xyz_direction_orig