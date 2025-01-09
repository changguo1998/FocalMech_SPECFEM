!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  subroutine compute_seismograms()

!   use constants, only: IOUT,OUTPUT_FILES
  ! seismograms
!   use specfem_par, only: NDIM,NGLOB_AB,nlength_seismogram,seismo_current,myrank, &
!     displ_record,veloc_record,accel_record
!   use specfem_par, only: NDIM,NGLOB_AB,NUsedGlob,seismo_current, &
!     used_globid,displ_record,veloc_record,accel_record
    ! use specfem_par, only: NDIM,NUsedGlob,seismo_current, &
    !   used_globid,veloc_record

  ! wavefields
!   use specfem_par_elastic, only: veloc

  implicit none

!   integer :: i,j,k
!   character(len=64) :: fname

!   do j = 1,NGLOB_AB
!     do i = 1,NDIM
!       displ_record(i,j,seismo_current) = displ(i,j)
!       veloc_record(i,j,seismo_current) = veloc(i,j)
!       accel_record(i,j,seismo_current) = accel(i,j)
!     enddo
!   enddo
!   do k = 1,NUsedGlob
!     j = used_globid(k)
!     do i = 1,NDIM
        ! displ_record(i,k,seismo_current) = displ(i,j)
        ! veloc_record(i,k,seismo_current) = veloc(i,j)
        ! accel_record(i,k,seismo_current) = accel(i,j)
!     enddo
!   enddo
!   write(fname, "('/proc',i6.6,'veloc_record.bin')") myrank
!   if (seismo_current == 1) then
!     open(unit=IOUT, file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//fname(1:len_trim(fname)), &
!     form='unformatted', action='write', access='stream', status='replace', iostat=ier)
!     write(IOUT) NDIM,NGLOB_AB,nlength_seismogram
!   else
!     open(unit=IOUT, file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//fname(1:len_trim(fname)), &
!     form='unformatted', action='write', access='stream', status='old', position='append', iostat=ier)
!   endif
!   if (ier /= 0) call exit_mpi(myrank,'Error opening file: '// &
!     OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//fname(1:len_trim(fname)))

!   write(IOUT) veloc
!   close(IOUT)
  end subroutine compute_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_seismograms_strain_adjoint()

  use constants, only: myrank,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use specfem_par, only: SIMULATION_TYPE,NGLOB_AB,ibool, &
    deltat,DT,t0,NSTEP,it, &
    seismo_current,seismo_offset,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    ispec_selected_source, &
    number_receiver_global,nrec_local, &
    Mxx,Myy,Mzz,Mxy,Mxz,Myz,tshift_src,hdur_Gaussian, &
    hprime_xx,hprime_yy,hprime_zz, &
    hxir_store,hetar_store,hgammar_store, &
    hpxir_store,hpetar_store,hpgammar_store, &
    ELASTIC_SIMULATION

  use specfem_par, only: GPU_MODE, Mesh_pointer

  ! strain "seismogram" and source derivatives (adjoint simulations)
  use specfem_par, only: seismograms_eps, &
    Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der,sloc_der

  ! wavefield
  use specfem_par_elastic, only: ispec_is_elastic,displ

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element

  integer :: irec_local,irec,idx
  integer :: iglob,ispec,i,j,k

  ! adjoint locals
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM):: eps_s
  real(kind=CUSTOM_REAL),dimension(NDIM):: eps_m_s
  real(kind=CUSTOM_REAL):: stf_deltat
  double precision :: stf
  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar
  double precision :: hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)

  double precision, external :: comp_source_time_function

  ! checks if anything to do
  if (SIMULATION_TYPE /= 2) return
  if (.not. ELASTIC_SIMULATION) return

  ! for strain seismograms
  ! strain seismograms are not implemented yet on GPU, thus transfers wavefield to CPU and computes it here
  ! transfers displacement to the CPU
  if (GPU_MODE) call transfer_displ_from_device(NDIM*NGLOB_AB, displ, Mesh_pointer)

  ! adjoint strain seismogram
  ! current index in seismogram_eps
  idx = seismo_offset + seismo_current

  ! checks bounds
  if (idx < 1 .or. idx > NSTEP/NTSTEP_BETWEEN_OUTPUT_SAMPLE) &
    call exit_mpi(myrank,'Error: seismograms_eps has wrong current index')

  ! loops over local receivers
  do irec_local = 1,nrec_local
    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! spectral element in which the receiver is located
    ! adjoint "receivers" are located at CMT source positions
    ! note: we take here xi_source,.. when FASTER_RECEIVERS_POINTS_ONLY is set
    ispec = ispec_selected_source(irec)

    ! additional calculations for pure adjoint simulations
    ! computes derivatives of source parameters

    ! elastic wave field
    if (ispec_is_elastic(ispec)) then
      ! stores elements displacement field
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            displ_element(:,i,j,k) = displ(:,iglob)
          enddo
        enddo
      enddo

      ! gets local receiver interpolators
      ! (1-D Lagrange interpolators)
      hxir(:) = hxir_store(:,irec_local)
      hetar(:) = hetar_store(:,irec_local)
      hgammar(:) = hgammar_store(:,irec_local)

      ! gets derivatives of local receiver interpolators
      hpxir(:) = hpxir_store(:,irec_local)
      hpetar(:) = hpetar_store(:,irec_local)
      hpgammar(:) = hpgammar_store(:,irec_local)

      ! computes the integrated derivatives of source parameters (M_jk and X_s)
      call compute_adj_source_frechet(ispec,displ_element, &
                                      Mxx(irec),Myy(irec),Mzz(irec), &
                                      Mxy(irec),Mxz(irec),Myz(irec), &
                                      eps_s,eps_m_s, &
                                      hxir,hetar,hgammar,hpxir,hpetar,hpgammar, &
                                      hprime_xx,hprime_yy,hprime_zz)

      ! stores strain value
      seismograms_eps(:,:,irec_local,idx) = eps_s(:,:)

      ! source time function value
      stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-tshift_src(irec),hdur_Gaussian(irec))

      stf_deltat = real(stf * deltat * NTSTEP_BETWEEN_OUTPUT_SAMPLE,kind=CUSTOM_REAL)

      ! integrated moment tensor derivatives
      Mxx_der(irec_local) = Mxx_der(irec_local) + eps_s(1,1) * stf_deltat
      Myy_der(irec_local) = Myy_der(irec_local) + eps_s(2,2) * stf_deltat
      Mzz_der(irec_local) = Mzz_der(irec_local) + eps_s(3,3) * stf_deltat
      Mxy_der(irec_local) = Mxy_der(irec_local) + 2 * eps_s(1,2) * stf_deltat
      Mxz_der(irec_local) = Mxz_der(irec_local) + 2 * eps_s(1,3) * stf_deltat
      Myz_der(irec_local) = Myz_der(irec_local) + 2 * eps_s(2,3) * stf_deltat

      ! source location derivative
      sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_s(:) * stf_deltat

    endif ! elastic

  enddo ! nrec_local

  end subroutine compute_seismograms_strain_adjoint

  subroutine compute_interpolated_partial_v_partial_x(veloc,NGLOB_AB, &
    ispec,NSPEC_AB,ibool, &
    hxir,hetar,hgammar, &
    dvxdx,dvydx,dvzdx,dvxdy,dvydy,dvzdy,dvxdz,dvydz,dvzdz)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ZERO

implicit none

double precision,intent(out) :: dvxdx,dvydx,dvzdx,dvxdy,dvydy,dvzdy,dvxdz,dvydz,dvzdz
integer,intent(in) :: ispec
integer,intent(in) :: NSPEC_AB,NGLOB_AB
real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(in) :: veloc

double precision,dimension(NGLLX),intent(in) :: hxir
double precision,dimension(NGLLY),intent(in) :: hetar
double precision,dimension(NGLLZ),intent(in) :: hgammar
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! local parameters
real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: dvxdxi,dvydxi,dvzdxi
double precision :: hlagrange
integer :: i,j,k


call compute_gradient_in_acoustic(ispec,veloc(1,:),dvxdxi)
call compute_gradient_in_acoustic(ispec,veloc(2,:),dvydxi)
call compute_gradient_in_acoustic(ispec,veloc(3,:),dvzdxi)

! perform the general interpolation using Lagrange polynomials
dvxdx = ZERO
dvydx = ZERO
dvzdx = ZERO
dvxdy = ZERO
dvydy = ZERO
dvzdy = ZERO
dvxdz = ZERO
dvydz = ZERO
dvzdz = ZERO

! interpolates seismograms at exact receiver locations
do k = 1,NGLLZ
do j = 1,NGLLY
do i = 1,NGLLX

hlagrange = hxir(i) * hetar(j) * hgammar(k)

dvxdx = dvxdx + dble(dvxdxi(1,i,j,k)) * hlagrange
dvydx = dvydx + dble(dvydxi(1,i,j,k)) * hlagrange
dvzdx = dvzdx + dble(dvzdxi(1,i,j,k)) * hlagrange
dvxdy = dvxdy + dble(dvxdxi(2,i,j,k)) * hlagrange
dvydy = dvydy + dble(dvydxi(2,i,j,k)) * hlagrange
dvzdy = dvzdy + dble(dvzdxi(2,i,j,k)) * hlagrange
dvxdz = dvxdz + dble(dvxdxi(3,i,j,k)) * hlagrange
dvydz = dvydz + dble(dvydxi(3,i,j,k)) * hlagrange
dvzdz = dvzdz + dble(dvzdxi(3,i,j,k)) * hlagrange

enddo
enddo
enddo


end subroutine compute_interpolated_partial_v_partial_x
