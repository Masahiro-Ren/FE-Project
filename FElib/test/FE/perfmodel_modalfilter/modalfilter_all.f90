#include "scalelib.h"

program perf_modalfilter
    use scale_precision
    use scale_element_line
    use scale_element_hexahedral
    use scale_element_modalfilter
    implicit none

    type(HexahedralElement) :: elem
    type(ModalFilter) :: filter

    integer, parameter :: Ne=342
    integer, parameter :: porder=11
    integer, parameter :: Np1D=porder+1
    integer, parameter :: Np=(porder+1)**3

    integer :: x, y, z, kp, ke, steps
    real(RP) :: DDENS(Np, Ne)
    real(RP) :: MOMX(Np, Ne)
    real(RP) :: MOMY(Np, Ne)
    real(RP) :: MOMZ(Np, Ne)
    real(RP) :: DRHOT(Np, Ne)
    real(RP) :: Gsqrt(Np, Ne)

    call elem%Init(porder, porder, .false.)

    call filter%Init(elem, 0D0, 0D0, 16, 0D0, 0D0, 16)

    !! Init
    !$omp parallel do private(x, y, z, kp)
    do ke = 1, Ne
        do z = 1, Np1D
        do y = 1, Np1D
        do x = 1, Np1D
                kp = x + (y - 1) * Np1D + (z - 1) * Np1D * Np1D
                DDENS(kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
                MOMX (kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
                MOMY (kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
                MOMZ (kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
                DRHOT(kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
                Gsqrt(kp,ke) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
        end do
        end do
        end do
    end do

    call fapp_start("modal_test", 1, 0)
    do steps = 1, 150
        call apply_modalfilter(Np, Ne, DDENS, MOMX, MOMY, MOMZ, DRHOT, filter, Gsqrt)
    end do
    call fapp_stop("modal_test", 1, 0)

    write(*,*) "Finished Peacefully."
        
   
contains
    !OCL SERIAL
    subroutine apply_modalfilter(Np, Ne, DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, filter, Gsqrt)
        implicit none

        integer, intent(in) :: Np, Ne
        real(RP), intent(inout) :: DDENS_(Np, Ne)
        real(RP), intent(inout) :: MOMX_(Np, Ne)
        real(RP), intent(inout) :: MOMY_(Np, Ne)
        real(RP), intent(inout) :: MOMZ_(Np, Ne)
        real(RP), intent(inout) :: DRHOT_(Np, Ne)
        class(ModalFilter), intent(in) :: filter
        real(RP), intent(in) :: Gsqrt(Np, Ne)


        integer :: ke, ii
        real(RP) :: tmp(Np, 5)
        real(RP) :: RGsqrt(Np)

        !$omp parallel do private(ii,tmp,RGsqrt)
        do ke = 1, Ne
            tmp(:,:) = 0.0_RP
            do ii=1, Np
                DDENS_(ii,ke) = Gsqrt(ii,ke) * DDENS_(ii,ke)
                MOMX_(ii,ke) = Gsqrt(ii,ke) * DDENS_(ii,ke)
                MOMY_(ii,ke) = Gsqrt(ii,ke) * DDENS_(ii,ke)
                MOMZ_(ii,ke) = Gsqrt(ii,ke) * DDENS_(ii,ke)
                DRHOT_(ii,ke) = Gsqrt(ii,ke) * DDENS_(ii,ke)
            end do

            call apply_filter_xyz_direction_split(filter%filterMat, filter%filterMat_tr, filter%filterMat_v_tr, DDENS_(:,ke), tmp(:,1))
            call apply_filter_xyz_direction_split(filter%filterMat, filter%filterMat_tr, filter%filterMat_v_tr, MOMX_ (:,ke), tmp(:,1))
            call apply_filter_xyz_direction_split(filter%filterMat, filter%filterMat_tr, filter%filterMat_v_tr, MOMY_ (:,ke), tmp(:,1))
            call apply_filter_xyz_direction_split(filter%filterMat, filter%filterMat_tr, filter%filterMat_v_tr, MOMZ_ (:,ke), tmp(:,1))
            call apply_filter_xyz_direction_split(filter%filterMat, filter%filterMat_tr, filter%filterMat_v_tr, DRHOT_(:,ke), tmp(:,1))

            RGsqrt(:) = 1.0_RP / Gsqrt(:,ke)
            DDENS_(:,ke) = tmp(:,1) * RGsqrt(:)
            MOMX_ (:,ke) = tmp(:,2) * RGsqrt(:)
            MOMY_ (:,ke) = tmp(:,3) * RGsqrt(:)
            MOMZ_ (:,ke) = tmp(:,4) * RGsqrt(:)
            DDENS_(:,ke) = tmp(:,5) * RGsqrt(:)
        end do

    end subroutine apply_modalfilter


    !OCL SERIAL
    subroutine apply_filter_xyz_direction_split(filterMat, filterMat_tr, filterMat_vtr, q_in, q_tmp )
        implicit none

        real(RP), intent(in) :: filterMat(12, 12)
        real(RP), intent(in) :: filterMat_tr(12, 12)
        real(RP), intent(in) :: filterMat_vtr(12, 12)
        real(RP), intent(inout) :: q_in(12, 12, 12)
        real(RP), intent(inout) :: q_tmp(12, 12, 12)
        
        real(RP) :: tmp1, tmp2, tmp3
        integer :: i, j, k

        !-- x direction
        do k=1, 12
        do j=1, 12
        do i=1, 12

        tmp1 = filterMat(i,1)  * q_in(1,j,k) + &
                filterMat(i,2)  * q_in(2,j,k) + &
                filterMat(i,3)  * q_in(3,j,k) + & 
                filterMat(i,4)  * q_in(4,j,k)
                
        tmp2 = filterMat(i,5)  * q_in(5,j,k) + & 
                filterMat(i,6)  * q_in(6,j,k) + & 
                filterMat(i,7)  * q_in(7,j,k) + & 
                filterMat(i,8)  * q_in(8,j,k) 

        tmp3 = filterMat(i,9)  * q_in(9,j,k) + & 
                filterMat(i,10) * q_in(10,j,k) + & 
                filterMat(i,11) * q_in(11,j,k) + & 
                filterMat(i,12) * q_in(12,j,k)  

        q_tmp(i,j,k) = tmp1 + tmp2 + tmp3

        end do
        end do
        end do

        !-- y direction
        do k=1, 12
        do j=1, 12
        do i=1, 12

        tmp1 = q_tmp(i,1,k)  * filterMat_tr(1,j) + &
                q_tmp(i,2,k)  * filterMat_tr(2,j) + &
                q_tmp(i,3,k)  * filterMat_tr(3,j) + &
                q_tmp(i,4,k)  * filterMat_tr(4,j)

        tmp2 = q_tmp(i,5,k)  * filterMat_tr(5,j) + &
                q_tmp(i,6,k)  * filterMat_tr(6,j) + &
                q_tmp(i,7,k)  * filterMat_tr(7,j) + &
                q_tmp(i,8,k)  * filterMat_tr(8,j)

        tmp3 = q_tmp(i,9,k)  * filterMat_tr(9,j) + &
                q_tmp(i,10,k) * filterMat_tr(10,j) + &
                q_tmp(i,11,k) * filterMat_tr(11,j) + &
                q_tmp(i,12,k) * filterMat_tr(12,j) 

        q_in(i,j,k) = tmp1 + tmp2 + tmp3

        end do
        end do
        end do

        !-- z direction
        do k=1, 12
        do j=1, 12
        do i=1, 12

        tmp1 = q_in(i,j,1)  * filterMat_vtr(1,k) + &
                q_in(i,j,2)  * filterMat_vtr(2,k) + &
                q_in(i,j,3)  * filterMat_vtr(3,k) + &
                q_in(i,j,4)  * filterMat_vtr(4,k)

        tmp2 = q_in(i,j,5)  * filterMat_vtr(5,k) + &
                q_in(i,j,6)  * filterMat_vtr(6,k) + &
                q_in(i,j,7)  * filterMat_vtr(7,k) + &
                q_in(i,j,8)  * filterMat_vtr(8,k)

        tmp3 = q_in(i,j,9)  * filterMat_vtr(9,k) + &
                q_in(i,j,10) * filterMat_vtr(10,k) + &
                q_in(i,j,11) * filterMat_vtr(11,k) + &
                q_in(i,j,12) * filterMat_vtr(12,k) 

        q_tmp(i,j,k) = tmp1 + tmp2 + tmp3

        end do
        end do
        end do
    end subroutine apply_filter_xyz_direction_split

    !OCL SERIAL
    subroutine apply_filter_xyz_direction_orig(filterMat, filterMat_tr, filterMat_vtr, q_in, q_tmp )
        implicit none

        real(RP), intent(in) :: filterMat(12, 12)
        real(RP), intent(in) :: filterMat_tr(12, 12)
        real(RP), intent(in) :: filterMat_vtr(12, 12)
        real(RP), intent(inout) :: q_in(12, 12, 12)
        real(RP), intent(inout) :: q_tmp(12, 12, 12)
        
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
end program perf_modalfilter
