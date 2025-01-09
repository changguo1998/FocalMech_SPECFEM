!---------------------------------------------------
! add by Guo Chang
! write out information to interpolate wavefield to receiver using other programs
!---------------------------------------------------

subroutine write_info_for_interpolation()
    use constants
    use specfem_par

    implicit none

    character(len=MAX_STRING_LEN) :: outputname
    integer :: fileid

    call synchronize_all()
    fileid = IOUT + myrank

    ! binary format
    if (USE_BINARY_FOR_SEISMOGRAMS) then
        if (myrank==0) then
            write(outputname, "('/globinterp.bin')")
            open(unit=fileid,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
                        form='unformatted',action='write',access='stream')
            write(fileid) NDIM,NGLOB_AB,NSPEC_AB,nlength_seismogram
            write(fileid) NGLLX,NGLLY,NGLLZ,nrec ! 10 int
            write(fileid) ibool ! NGLLX x NGLLY x NGLLZ x NSPEC
            write(fileid) ispec_selected_rec
            write(fileid) xixstore,xiystore,xizstore! NGLLX x NGLLY x NGLLZ x NSPEC
            write(fileid) etaxstore,etaystore,etazstore! NGLLX x NGLLY x NGLLZ x NSPEC
            write(fileid) gammaxstore,gammaystore,gammazstore! NGLLX x NGLLY x NGLLZ x NSPEC
            write(fileid) network_name,station_name
            write(fileid) hprime_xx,hprime_yy,hprime_zz ! GLL x GLL
            write(fileid) nu_rec
            close(fileid)
        endif
        if (nrec > 2) then
            write(outputname, "('/proc',i6.6,'localinterp.bin')") myrank
            open(unit=fileid,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
                        form='unformatted',action='write',access='stream')
            write(fileid) NDIM,NGLLX,nlength_seismogram,nrec_local,NUsedGlob
            write(fileid) used_globid
            write(fileid) number_receiver_global ! LREC
            write(fileid) hxir_store,hetar_store,hgammar_store  !NGLL x LREC
            close(fileid)
            write(outputname, "('/proc',i6.6,'specinfo.bin')") myrank
            open(unit=fileid,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
                        form='unformatted',action='write',access='stream')
            write(fileid) NUsedGlob
            write(fileid) used_globid
            close(fileid)
        endif
    endif


end subroutine write_info_for_interpolation

! subroutine write_info_for_interpolation()

!     use constants
!     use specfem_par

!     implicit none

!     character(len=MAX_STRING_LEN) :: outputname
!     integer :: fileid,ispec,ispecir,iglob,lspec,lglob,irl,irg,i,j,k,nglobused,nspecused
!     logical,dimension(NGLOB_AB) :: glob_is_used
!     integer,dimension(NGLOB_AB) :: used_glob_id
!     logical,dimension(NSPEC_AB) :: spec_is_used
!     integer,dimension(NSPEC_AB) :: used_spec_id

!     integer,dimension(:),allocatable :: lspec_select_recl
!     integer,dimension(:,:,:,:),allocatable :: ibool_l
!     real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: xixl,xiyl,xizl, &
!         etaxl,etayl,etazl,gammaxl,gammayl,gammazl
!     character(len=MAX_LENGTH_STATION_NAME),dimension(:),allocatable :: stationl
!     character(len=MAX_LENGTH_NETWORK_NAME),dimension(:),allocatable :: networkl
!     double precision,dimension(:,:,:),allocatable :: nu_recl
!     real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: velocl
!     call synchronize_all()
!     fileid = IOUT + myrank

!     do iglob = 1,NGLOB_AB
!         glob_is_used(iglob) = .false.
!     enddo
!     do ispec = 1,NSPEC_AB
!         spec_is_used(ispec) = .false.
!     enddo

!     do irl = 1,nrec_local
!         irg = number_receiver_global(irl)
!         ispec = ispec_selected_rec(irg)
!         spec_is_used(ispec) = .true.
!         do k = 1,NGLLZ
!             do j = 1,NGLLY
!                 do i = 1,NGLLX
!                     iglob = ibool(i,j,k,ispec)
!                     glob_is_used(iglob) = .true.
!                 enddo
!             enddo
!         enddo
!     enddo
!     nspecused = 0
!     do ispec = 1,NSPEC_AB
!         if (spec_is_used(ispec)) then
!             nspecused = nspecused + 1
!             used_spec_id(nspecused) = ispec
!         endif
!     enddo
!     nglobused = 0
!     do iglob = 1,NGLOB_AB
!         if (glob_is_used(iglob)) then
!             nglobused = nglobused + 1
!             used_glob_id(nglobused) = iglob
!         endif
!     enddo

!     allocate(lspec_select_recl(nrec_local))
!     allocate(ibool_l(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(xixl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(xiyl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(xizl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(etaxl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(etayl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(etazl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(gammaxl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(gammayl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(gammazl(NGLLX,NGLLY,NGLLZ,nspecused))
!     allocate(networkl(nrec_local))
!     allocate(stationl(nrec_local))
!     allocate(nu_recl(NDIM,NDIM,nrec_local))
!     allocate(velocl(NDIM,nglobused,nlength_seismogram))

!     do irl = 1,nrec_local
!         irg = number_receiver_global(irl)
!         ispec = ispec_selected_rec(irg)
!         do lspec = 1,nspecused
!             if (used_spec_id(lspec) == ispec) then
!                 lspec_select_recl(irl) = lspec
!             endif
!         enddo
!         networkl(irl) = network_name(irg)
!         stationl(irl) = station_name(irg)
!         do j = 1,NDIM
!             do i = 1,NDIM
!                 nu_recl(i,j,irl) = nu_rec(i,j,irg)
!             enddo
!         enddo
!     enddo

!     do lspec = 1,nspecused
!         ispec = used_spec_id(lspec)
!         ispecir = irregular_element_number(ispec)
!         do k = 1,NGLLZ
!         do j = 1,NGLLY
!         do i = 1,NGLLX
!             iglob = ibool(i,j,k,ispec)
!             do lglob = 1,nglobused
!                 if (used_glob_id(lglob) == iglob) then
!                     ibool_l(i,j,k,lspec) = lglob
!                 endif
!             enddo
!             if (ispecir > 0) then
!                 xixl(i,j,k,lspec) = xixstore(i,j,k,ispecir)
!                 xiyl(i,j,k,lspec) = xiystore(i,j,k,ispecir)
!                 xizl(i,j,k,lspec) = xizstore(i,j,k,ispecir)
!                 etaxl(i,j,k,lspec) = etaxstore(i,j,k,ispecir)
!                 etayl(i,j,k,lspec) = etaystore(i,j,k,ispecir)
!                 etazl(i,j,k,lspec) = etazstore(i,j,k,ispecir)
!                 gammaxl(i,j,k,lspec) = gammaxstore(i,j,k,ispecir)
!                 gammayl(i,j,k,lspec) = gammaystore(i,j,k,ispecir)
!                 gammazl(i,j,k,lspec) = gammazstore(i,j,k,ispecir)
!             else
!                 xixl(i,j,k,lspec) = xix_regular
!                 xiyl(i,j,k,lspec) = xix_regular
!                 xizl(i,j,k,lspec) = xix_regular
!                 etaxl(i,j,k,lspec) = xix_regular
!                 etayl(i,j,k,lspec) = xix_regular
!                 etazl(i,j,k,lspec) = xix_regular
!                 gammaxl(i,j,k,lspec) = xix_regular
!                 gammayl(i,j,k,lspec) = xix_regular
!                 gammazl(i,j,k,lspec) = xix_regular
!             endif
!         enddo
!         enddo
!         enddo
!     enddo
!     do lglob = 1,nglobused
!         iglob = used_glob_id(lglob)
!         do j = 1,nlength_seismogram
!             do i = 1,NDIM
!                 velocl(i,lglob,j) = veloc_record(i,iglob,j)
!             enddo
!         enddo
!     enddo

!     ! binary format
!     if (USE_BINARY_FOR_SEISMOGRAMS) then
!         write(outputname, "('/proc',i6.6,'interp.bin')") myrank
!         open(unit=fileid,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
!                     form='unformatted',action='write',access='stream')

!         write(fileid) myrank,NDIM,NGLOB_AB,NSPEC_AB,nglobused,nspecused,nlength_seismogram,NGLLX,NGLLY,NGLLZ,nrec,nrec_local ! 9 int
!         write(fileid) number_receiver_global ! LREC
!         write(fileid) lspec_select_recl     ! LREC
!         write(fileid) ibool_l               ! NGLLX x NGLLY x NGLLZ x LSPEC
!         write(fileid) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl ! NGLLX x NGLLY x NGLLZ x LSPEC
!         write(fileid) hprime_xx,hprime_yy,hprime_zz ! GLL x GLL
!         write(fileid) networkl,stationl     ! LREC
!         write(fileid) nu_recl               ! NDIM x NDIM x LREC
!         write(fileid) hxir_store,hetar_store,hgammar_store  !NGLL x LREC
!         write(fileid) velocl                ! NDIM x LGLOB x NT

!         close(fileid)

!     else
!         write(outputname, "('/proc',i6.6,'interp.ascii')") myrank
!         open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
!                     form='formatted',action='write')
!         write(fileid,*) myrank
!         write(fileid,*) NDIM
!         write(fileid,*) nglobused
!         write(fileid,*) nspecused
!         write(fileid,*) nlength_seismogram
!         write(fileid,*) NGLLX
!         write(fileid,*) NGLLY
!         write(fileid,*) NGLLZ
!         write(fileid,*) nrec_local
!         do irl = 1,nrec_local
!             write(fileid,*) lspec_select_recl(irl)     ! LREC
!         enddo
!         call print_4d_int(fileid,ibool_l,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,xixl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,xiyl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,xizl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,etaxl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,etayl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,etazl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,gammaxl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,gammayl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_4d(fileid,gammazl,NGLLX,NGLLY,NGLLZ,nspecused)
!         call print_2d(fileid,hprime_xx,NGLLX,NGLLX)
!         call print_2d(fileid,hprime_yy,NGLLY,NGLLY)
!         call print_2d(fileid,hprime_zz,NGLLZ,NGLLZ)
!         do irl = 1,nrec_local
!             write(fileid,*) networkl(irl)
!         enddo
!         do irl = 1,nrec_local
!             write(fileid,*) stationl(irl)
!         enddo
!         call print_3d_double(fileid,nu_recl,NDIM,NDIM,nrec_local)
!         call print_2d(fileid,hxir_store,NGLLX,nrec_local)
!         call print_2d(fileid,hetar_store,NGLLY,nrec_local)
!         call print_2d(fileid,hgammar_store,NGLLZ,nrec_local)
!         call print_3d(fileid,velocl,NDIM,nglobused,nlength_seismogram)
!         close(IOUT)
!     endif

!     deallocate(lspec_select_recl)
!     deallocate(ibool_l)
!     deallocate(xixl)
!     deallocate(xiyl)
!     deallocate(xizl)
!     deallocate(etaxl)
!     deallocate(etayl)
!     deallocate(etazl)
!     deallocate(gammaxl)
!     deallocate(gammayl)
!     deallocate(gammazl)
!     deallocate(networkl)
!     deallocate(stationl)
!     deallocate(nu_recl)
!     deallocate(velocl)

! end subroutine


subroutine print_1d(io,v,n)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n
    real(kind=CUSTOM_REAL),dimension(n),intent(in) :: v

    integer :: i

    do i = 1,n
        write(io,*) v(i)
    enddo
    write(io,*) '='
end subroutine print_1d

subroutine print_1d_int(io,v,n)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n
    integer,dimension(n),intent(in) :: v

    integer :: i

    do i = 1,n
        write(io,*) v(i)
    enddo
    write(io,*) '='
end subroutine print_1d_int

subroutine print_2d(io,v,n1,n2)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n1,n2
    real(kind=CUSTOM_REAL),dimension(n1,n2),intent(in) :: v

    integer :: i,j

    do j=1,n2
        do i=1,n1
            write(io,*) v(i,j)
        enddo
    enddo
    write(io,*) '='
end subroutine print_2d

subroutine print_3d(io,v,n1,n2,n3)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n1,n2,n3
    real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(in) :: v

    integer :: i,j,k

    do k=1,n3
    do j=1,n2
    do i=1,n1
        write(io,*) v(i,j,k)
    enddo
    enddo
    enddo
    write(io,*) '='
end subroutine print_3d

subroutine print_3d_double(io,v,n1,n2,n3)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n1,n2,n3
    double precision,dimension(n1,n2,n3),intent(in) :: v

    integer :: i,j,k

    do k=1,n3
    do j=1,n2
    do i=1,n1
        write(io,*) v(i,j,k)
    enddo
    enddo
    enddo
    write(io,*) '='
end subroutine print_3d_double

subroutine print_4d(io,v,n1,n2,n3,n4)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n1,n2,n3,n4
    real(kind=CUSTOM_REAL),dimension(n1,n2,n3,n4),intent(in) :: v

    integer :: i,j,k,l

    do l=1,n4
    do k=1,n3
        do j=1,n2
        do i=1,n1
            write(io,*) v(i,j,k,l)
        enddo
        enddo
    enddo
    enddo
    write(io,*) '='
end subroutine print_4d

subroutine print_4d_int(io,v,n1,n2,n3,n4)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: io,n1,n2,n3,n4
    integer,dimension(n1,n2,n3,n4),intent(in) :: v

    integer :: i,j,k,l

    do l=1,n4
    do k=1,n3
        do j=1,n2
        do i=1,n1
            write(io,*) v(i,j,k,l)
        enddo
        enddo
    enddo
    enddo
    write(io,*) '='
end subroutine print_4d_int
