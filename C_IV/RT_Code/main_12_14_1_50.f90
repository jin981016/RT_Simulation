program radiative_transfer
use random
use mpi
use cons
!use obs_mod
use observer_mod
use escape_obs_mod
use tracing_tau_mod
use gen_photon_mod
use grid_mod
use scattering_mod
use peeling_off_mod
use char_mod
use atom_mod
use clumps_mod
use clumpy_RT_mod
use dust_mod
implicit none

integer, parameter :: nphoton_emit = 1.0e7
integer :: nphoton_flat
!	MAIN
integer, parameter :: nN_atom = 32, nv_exp = 17, nv_ran = 17, ntau_d = 17
real(kind=rkd) :: N_atom(nN_atom), v_exp(nv_exp), v_ran(nv_ran), tau_d(ntau_d+1)
real(kind=rkd) :: N_atom_min, N_atom_max, dN_atom
real(kind=rkd) :: tau_d_min, tau_d_max, dtau_d
integer :: iN_atom, iv_exp, iv_ran, itau_d
!	clumps parameter
integer, parameter :: nf_c = 11, nr_cl = 1, nv_cl = 5
real(kind=rkd) :: f_c(nf_c), r_cl(nr_cl), v_cl(nv_cl)
integer :: if_c, ir_cl, iv_cl
real(kind=rkd) :: atom_frac
!
integer, parameter :: nv_emit = 17
real(kind=rkd), dimension(nv_emit) :: v_emit
integer :: iv_emit
!
integer, parameter :: nN_atom_cl = 3, nv_exp_cl = 1, nv_ran_cl = 3
real(kind=rkd) :: N_atom_cl(nN_atom_cl), v_exp_cl(nv_exp_cl), v_ran_cl(nv_ran_cl)
integer :: iN_atom_cl, iv_exp_cl, iv_ran_cl
integer, parameter :: nR_out = 5
real(kind=rkd) :: R_out(nR_out)
integer :: iR_out
integer :: iphoton
integer :: NS
logical :: temp_phi
!	MPI
integer :: nproc, rank, irank, ip
integer :: err, num_photon = 0
integer :: ans , tag
integer :: anstype, sender 
integer :: status(MPI_STATUS_SIZE), request_count
real(kind=rkd) :: time1, time2
real(kind=rkd) :: tau_s 
type(photon_type) :: photon
type(photon_type) :: photon_i

CALL MPI_INIT(err)

call MPI_COMM_SIZE(MPI_COMM_WORLD,mpar%nproc,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%rank,ierr)
call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpar%h_comm, ierr)

! slited by each node
call MPI_COMM_RANK(mpar%h_comm, mpar%h_rank,ierr)
call MPI_COMM_SIZE(mpar%h_comm, mpar%h_nproc, ierr)


! slited by same rank in each node
call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%rank, mpar%SUB_COMM,ierr)
call MPI_COMM_SIZE(mpar%SUB_COMM, mpar%sub_nproc,ierr)





call set_atoms()
atom = C_IV
 
!v_ran(1) = 50.d5
!v_ran(2) = 70.d5
!v_ran(3) = 80.d5
!v_ran(4) = 90.d5
!v_ran(1) = 100.d5

 
!v_ran(1) = sqrt(2.d0*k*1.e5/atom%mass)
v_ran(1) = 50.d5
v_ran(2) = 150.d5
v_ran(3) = 250.d5
v_ran(4) = 350.d5
v_ran(5) = 450.d5
v_ran(6) = 550.d5
v_ran(7) = 650.d5
v_ran(8) = 750.d5
v_ran(9) = 850.d5
v_ran(10) = 950.d5


iv_ran = 1

tau_D(1) = 0.d0


v_exp(1)  = 50.d5	! cm/s
v_exp(2)  = 150.d5	! cm/s
v_exp(3)  = 250.d5	! cm/s
v_exp(4)  = 350.d5
v_exp(5)  = 450.d5	! cm/s
v_exp(6)  = 550.d5	! cm/s
v_exp(7)  = 650.d5 
v_exp(8)  = 750.d5      ! cm/s
v_exp(9)  = 850.d5      ! cm/s
v_exp(10)  = 950.d5
v_exp(11)  = 1000.d5

v_emit(1)  = 50.d5       ! cm/s
v_emit(2)  = 150.d5     ! cm/s
v_emit(3)  = 250.d5     ! cm/s
v_emit(4)  = 350.d5
v_emit(5)  = 450.d5     ! cm/s
v_emit(6)  = 550.d5     ! cm/
v_emit(7)  = 650.d5
v_emit(8)  = 750.d5       ! cm/s
v_emit(9)  = 850.d5     ! cm/s
v_emit(10)  = 950.d5     ! cm/s
v_emit(11)  = 1000.d5     ! cm/s


N_atom(1) = 1.3d12	! cm^-2
N_atom(2) = 2.d12	! cm^-2
N_atom(3) = 3.2d12	! cm^-2
N_atom(4) = 5.d12	! cm^-2
N_atom(5) = 7.9d12	! cm^-2


N_atom(6) = 1.3d13      ! cm^-2
N_atom(7) = 2.d13      ! cm^-2
N_atom(8) = 3.2d13      ! cm^-2
N_atom(9) = 5.d13       ! cm^-2
N_atom(10) = 7.9d13      ! cm^-2


N_atom(11) = 1.3d14	! cm^-2
N_atom(12) = 2.d14	! cm^-2
N_atom(13) = 3.2d14	! cm^-2
N_atom(14) = 5.d14	! cm^-2
N_atom(15) = 7.9d14	! cm^-2

N_atom(16) = 1.3d15	! cm^-2
N_atom(17) = 2.d15	! cm^-2
N_atom(18) = 3.2d15	! cm^-2
N_atom(19) = 5.d15	! cm^-2
N_atom(20) = 7.9d15	! cm^-2

N_atom(21) = 1.3d16	! cm^-2
N_atom(22) = 2.d16	! cm^-2
N_atom(23) = 3.2d16	! cm^-2
N_atom(24) = 5.d16	! cm^-2
N_atom(25) = 7.9d16	! cm^-2

N_atom(26) = 1.3d17	! cm^-2
N_atom(27) = 2.0d17
N_atom(28) = 3.2d17	! cm^-2
N_atom(29) = 5.d17
N_atom(30) = 7.9d17



itau_d = 1


do iv_emit = 1,10
do iv_ran = 1,10
do iv_exp =1,5
do iN_atom = 1,15




call set_escape_observer()
call set_dust('dust_data/MW_C_IV.dat')
	write(fn_model,100) 'data_real/N_atom',N_atom(iN_atom), &
					'_Vexp', v_exp(iv_exp)/1e5, &
					'_Vemit', v_emit(iv_emit)/1e5, &
					'_tauD', tau_d(itau_d), &
					'_Vran', v_ran(iv_ran)/1e5
100 	format(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)
        call trim_char(fn_model)


	
if(mpar%rank .eq. master)	print*,fn_model



if(mpar%rank .eq. master) then
print*,'Start'
endif


time1 = MPI_WTIME()

atom = C_IV
!call set_observer()
call init_random_seed()
call MPI_BARRIER(MPI_COMM_WORLD,err)

!call set_grid_empty
call set_grid_sphere(N_atom(iN_atom), v_ran(iv_ran), v_exp(iv_exp), tau_d(itau_d))
if(mpar%rank .eq. master) then
print*,'set_grid_done'
endif
call set_clumps_no()
!call set_clumps(f_c(if_c), N_atom(iN_atom), v_exp(iv_exp), r_cl(ir_cl), v_ran(iv_ran), tau_d(itau_d), v_cl(iv_cl) )

if(mpar%rank .eq. master) then
print*,'set_clumps_done'
endif


nphoton = nphoton_emit
num_photon = 0


if(mpar%rank .eq. master) then
print*,'Done'
time2 = MPI_WTIME()
print*,'TIME : ',time2 - time1, ' sec'
endif

call MPI_BARRIER(MPI_COMM_WORLD,err)

if(mpar%rank .eq. master) then	! Master

	do irank = 1, min(nphoton, mpar%nproc-1) 
	ip = irank
	tag = 1
	call MPI_SEND(ip,1,MPI_INTEGER8, irank, tag,MPI_COMM_WORLD,err)
	num_photon = num_photon + 1
	enddo

	do irank = 1, nphoton
	call MPI_RECV(ans,1,MPI_INTEGER8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,err)
	sender  = status(MPI_SOURCE)
	anstype = status(MPI_TAG)

	if(num_photon .lt. nphoton) then

	if(mod(num_photon,nphoton/2) .eq. 0) then
	time2 = MPI_WTIME()
	print*,num_photon*100/nphoton,'%',sender,time2-time1
	endif


	tag = 1
	num_photon = num_photon + 1
	ip = num_photon
	call MPI_SEND(ip,1,MPI_INTEGER8, sender, tag, MPI_COMM_WORLD,err)

	else

	tag = 0
	call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,sender,tag,MPI_COMM_WORLD,err)

	endif

	enddo

else
call init_random_seed()
	time_i = MPI_WTIME()
	do 607
        call MPI_RECV(ip,1,MPI_INTEGER8,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,err)
	photon%ip = ip
!	print*,rank,'slave'
        if (status(MPI_TAG) .eq.  0) exit
	ans = 1

	call gen_photon_Gaussian(photon,v_emit(iv_emit))
	call peeling_off_direct_metal(photon)

		do 608

                tau_s = -dlog(rand_number())
                call tracing_tau_clumpy(photon,tau_s)

                        if(photon%esc .eqv. .true.) then
                        exit
                        else
!	Scattering
				if(photon%clump .eqv. .false.) then

	                        call scattering(photon)

				else if(photon%clump .eqv. .true.) then
				call RT_in_clump(photon)

	                        if(photon%esc .eqv. .true.) exit

				endif
                        endif

608		continue


call collect_escape_photon(photon)

tag = 1
call MPI_SEND(ans,1,MPI_INTEGER8,master,tag,MPI_COMM_WORLD,err)
607	continue

endif	!	Slave


if(mpar%rank .eq. master) then
	time2 = MPI_WTIME()
	print*,'	 100 %',sender,time2-time1
endif


call MPI_BARRIER(MPI_COMM_WORLD,err)
call reduce_escape_observer()
call MPI_BARRIER(MPI_COMM_WORLD,err)
call write_escape_observer()
call MPI_BARRIER(MPI_COMM_WORLD,err)

call clear_grid()
call clear_escape_observer()
call clear_clumps()
call clear_dust()

call MPI_BARRIER(MPI_COMM_WORLD,err)

time2 = MPI_WTIME()

if(mpar%rank .eq. master) then
print*,'Time :',time2-time1
endif


enddo
enddo
enddo
enddo



call MPI_BARRIER(MPI_COMM_WORLD,err)
CALL MPI_FINALIZE(err)


end program radiative_transfer
