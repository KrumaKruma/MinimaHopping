
! program testf2p
!   implicit none
!   integer, parameter :: nat = 16
!   integer, parameter:: natx_sphere=40! = 35  ! maximum number of atoms in the sphere
!   integer, parameter :: ns =1, np=1 !  = 1 = 1
!   integer :: num_cat=1 !Number of atom enviroments per chemical element (e.g. two different enviroments for C atoms)
!   integer :: nex_cutoff=2
!   character(len=2):: atomnames(nat)
!   real*8, dimension(3, nat) :: rxyz !atomic input positions (unchange on output)
!   real*8, dimension(3, 3) :: alat ! three column vectors specifying the cell vectors
!   real*8 :: width_cutoff =4.0d0
!   integer :: lengthfp = 160
!   integer :: nat_sphere_current_max
!
!   real*8:: penalty
!   real*8, dimension(3, nat) :: dpenaldr ! derivative of penalty function with respect to all atomic positions
!   real*8 :: dpenaldalat(3,3)
!
!   character(len=150) :: PosFile
!
!   real*8 :: Bohr_Ang=0.529177d0  ! Conversion factor Bohr to angstroem
!
!   ! if (file_format_ascii_logical) then !.ascii case
!   PosFile = 'poscur.ascii'
!   ! else !.xyz type
!   !   PosFile = 'poscur.xyz'
!   ! end if
!
!   call read_ascii_extern(PosFile, rxyz, alat, nat, atomnames)
!
!   ! rxyz = rxyz / Bohr_Ang
!   ! alat = alat / Bohr_Ang
!
!
!   call SymmetryPenalty_v3(nat, alat, rxyz, atomnames, natx_sphere, ns, np, width_cutoff, num_cat, nex_cutoff, &
!                              lengthfp, nat_sphere_current_max, dpenaldr, penalty, dpenaldalat)
!
!   write(*,*) penalty, 'penalty'
!   write(*,*) 'dpenaldr'
!   write(*,*) dpenaldr
!   write(*,*) 'dpenaldalat'
!   write(*,*) dpenaldalat
!
! end program testf2p
!
!
! subroutine read_ascii_extern(filename, rxyz, alat, nat, atomnames)
!   implicit none
!   !! reads an ascii file with the specified filename units are assumed to be angstroem.
!   !! units are convertet to hartree units before they are returned
!   character(len=*), intent(in) :: filename
!   !! filename of the file which is created
!   integer, intent(in) :: nat
!   !! number of atoms
!   real*8, dimension(3, nat) :: rxyz
!   !! atom positions
!   real*8, dimension(3, 3) :: alat
!   !! lattice vectors
!   character(len=2) :: atomnames(nat)
!   !! String containing the chemical symbol of each atom.
!   integer :: i, io, ios
!   character(len=150) :: all_line
!   real*8 :: alat_temp(2, 3)
!   real*8 :: Bohr_Ang = 0.529177249d0
!
!   open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status="old")
!   if (ios /= 0) then
!     print *, "error opening ascii file: "//filename
!     stop
!   end if
!   read (io, *, iostat=ios) all_line
!   if (ios /= 0) then
!     print *, trim(adjustl(filename)), ios
!     stop "error reading file "
!   end if
!   read (io, *, iostat=ios) alat_temp(1, 1), alat_temp(1, 2), alat_temp(1, 3)
!   read (io, *, iostat=ios) alat_temp(2, 1), alat_temp(2, 2), alat_temp(2, 3)
!   if (ios /= 0) stop "reading lattice vectors"
!   alat = 0.0
!   alat(1, 1) = alat_temp(1, 1)
!   alat(1, 2) = alat_temp(1, 2)
!   alat(2, 2) = alat_temp(1, 3)
!   alat(1, 3) = alat_temp(2, 1)
!   alat(2, 3) = alat_temp(2, 2)
!   alat(3, 3) = alat_temp(2, 3)
!   i = 1
!   do while (i <= nat)
!     read (io, *, iostat=ios) rxyz(1, i), rxyz(2, i), rxyz(3, i), atomnames(i)
!     if (ios > 0) then
!       print *, rxyz(:, i), atomnames(i)
!       cycle
!     end if
!     if (ios < 0) then
!       print *, "end of file in read ascii, file: "//filename
!       stop
!     end if
!     i = i + 1
!   end do
!   close (io)
!   alat = alat/Bohr_Ang
!   rxyz = rxyz/Bohr_Ang
! end subroutine read_ascii_extern


! module symmetry_penalty
!   IMPLICIT NONE
! CONTAINS

! subroutine atomnames_test(nat, atom_element_number_array)
! implicit none
! integer, intent(in) :: nat
! integer, dimension(nat), intent(in) :: atom_element_number_array
! character(len=2) :: atomnames(nat)
! integer :: iat
!
! print*, nat ,'nat'
! print*, 'atom_element_number_array'
! print*, atom_element_number_array
!
! do iat = 1, nat
!   call element_number_2_atomname(atom_element_number_array(iat), atomnames(iat))
! end do
!
! ! print*, atomnames
!
! end subroutine atomnames_test


subroutine sympen_f2py_interface(nat, alat_in, rxyz_in, atom_element_number_array, natx_sphere, ns, np, width_cutoff, &
                           num_cat, nex_cutoff, lengthfp, nat_sphere_current_max, dpenaldr, penalty, dpenaldalat)

  implicit none
  integer, intent(in) :: nat
  real*8, dimension(3, 3), intent(in) :: alat_in ! three column vectors specifying the cell vectors
  real*8, dimension(3, 3) :: alat ! three column vectors specifying the cell vectors, intern
  real*8, dimension(3, nat), intent(in) :: rxyz_in !atomic input positions (unchange on output)
  real*8, dimension(3, nat) :: rxyz !atomic input positions, intern
  integer, dimension(nat), intent(in) :: atom_element_number_array
  integer, intent(in) :: natx_sphere! = 35  ! maximum number of atoms in the sphere
  integer, intent(in) :: ns, np !  = 1 = 1
  real*8, intent(in) :: width_cutoff
  integer, intent(in) :: num_cat !Number of atom enviroments per chemical element (e.g. two different enviroments for C atoms)
  integer, intent(in) :: nex_cutoff
  integer, intent(in) :: lengthfp
  integer, intent(out) :: nat_sphere_current_max
  integer :: num_dlamda !Number of calculated derivates of Evalues from DimensionalityMatrix starting from largest Evalue

  real*8, intent(out) :: penalty
  real*8, dimension(3, nat), intent(out) :: dpenaldr ! derivative of penalty function with respect to all atomic positions
  real*8, intent(out) :: dpenaldalat(3,3)

  character(len=2) :: atomnames(nat)
  integer :: iat
  real*8 :: Bohr_Ang = 0.529177249d0

  !Coversion of units done in python part
  alat = alat_in !/Bohr_Ang
  rxyz = rxyz_in !/Bohr_Ang

  ! print*, nat ,'nat'
  ! print*, 'atom_element_number_array'
  ! print*, atom_element_number_array

  do iat = 1, nat
    call element_number_2_atomname(atom_element_number_array(iat), atomnames(iat))
  end do

  ! print*, nex_cutoff, "nex_cutoff"

  call SymmetryPenalty_v3(nat, alat, rxyz, atomnames, natx_sphere, ns, np, width_cutoff, num_cat, nex_cutoff, &
                            lengthfp, nat_sphere_current_max, dpenaldr, penalty, dpenaldalat)

  ! print*, "penalty"
  ! print*, penalty


end subroutine sympen_f2py_interface



subroutine SymmetryPenalty_v3(nat, alat, rxyz, atomnames, natx_sphere, ns, np, width_cutoff, num_cat, nex_cutoff, &
                           lengthfp, nat_sphere_current_max, dpenaldr, penalty, dpenaldalat) !check_value, check_derdr, check_deralat
!  Calculates penalty from DimensionalityMatrix which gets calulated using the atomic fingerprints from xyz2devaldr
!  Also calculates the Derivate dpenaldr

!SymmetryPenalty needs rxyz in Bohr
use omp_lib
implicit none
integer, intent(in) :: nat
integer, intent(in) :: natx_sphere! = 35  ! maximum number of atoms in the sphere
integer, intent(in) :: ns, np !  = 1 = 1
integer, intent(in) :: num_cat !Number of atom enviroments per chemical element (e.g. two different enviroments for C atoms)
integer, intent(in) :: nex_cutoff
integer :: num_dlamda !Number of calculated derivates of Evalues from DimensionalityMatrix
!                                            starting from largest Evalue
! double precision, intent(out) :: check_value
! double precision, intent(out) :: check_derdr(3,nat)
! double precision, intent(out) :: check_deralat(3,3)
!! array containg the chemical symbols of each atom
character(len=2), intent(in) :: atomnames(nat)
real*8, dimension(3, nat), intent(in) :: rxyz !atomic input positions (unchange on output)
real*8, dimension(3, 3), intent(in) :: alat ! three column vectors specifying the cell vectors
real*8, intent(in) :: width_cutoff
integer, intent(in) :: lengthfp
integer, intent(out) :: nat_sphere_current_max
real*8, intent(out) :: penalty
real*8, dimension(3, nat), intent(out) :: dpenaldr ! derivative of penalty function with respect to all atomic positions
real*8, intent(out) :: dpenaldalat(3,3)
real*8, dimension(nat) :: rcov ! covalent radii neede for fingerprints
real*8, allocatable ::  rxyz_sphere(:, :), rcov_sphere(:) ! positions and covalent radii of atoms in sphere
real*8, allocatable ::  amplitude(:), deramplitude(:) ! amplitude modulationf function for Gaussians and its derivative
real*8, allocatable ::  fp_loc(:) ! eigenvlues of the overlap matric which give the finperprint vector for the atomic environment
real*8, allocatable ::  fp(:, :)   ! All the fingerprints for all the atoms
real*8, allocatable ::  dfpdr(:, :, :) ! dfpdr(ixyz,iat_sp,l)  derivative of l-th eigenvalues with respect
!                                       to the component ixyz  of the position of atom iat_sp in the sphere
real*8, allocatable :: dfp_marco(:,:,:,:)
!double precision ::  dFPdalatAll(nat, 3, 3, natx_sphere*(ns + np*3)) !dFPdrAll(iatsphere,ixyz,jlat,l)  derivatice of the l-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of lattice vector jlat
real*8, allocatable :: dFPdalatAll(:,:,:,:) !dFPdrAll(iatsphere,ixyz,jlat,l)  derivatice of the l-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of lattice vector jlat
real*8, allocatable ::  dfpdr0(:, :, :) ! Work array (same as dfpdr) but for a number of atoms in the sphere that is less or
!                                           equal to natx_sphere
!  dimension dlamdadr(3,nat,nat) ! dlamdadr(ixyz,iat,lat) is the derivative of Eval lat from DimensionalityMatrix with respect
!                                                 to all atomic positions
real*8, allocatable ::  dFPdrAll(:, :, :, :) !dFPdrAll(iat,ixyz,jat,iatsphere)  derivatice of the jat-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsohere, with respect to the component ixyz of atom iat
real*8, allocatable ::  DimMatrix(:, :) ! dimensionality matrix nat*nat
real*8, allocatable ::  evalsDimM(:)   ! eigenvalues of dimensionality matrix
real*8, allocatable ::  dDimMdr(:, :, :, :)   !dDimMdr(iat,ixyz,lat,jat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
real*8, allocatable ::  WORK(:)
real*8, allocatable ::  EigenvecDimM(:, :) ! can be eliminated
real*8, allocatable ::  tmpResdL(:, :, :) !  work array
integer, allocatable ::  indat(:) !contains  the indices of the atoms in the sphere
integer, allocatable ::  nat_sphere_array(:) ! needed to calculate maximum number of atoms in sphere
real*8 :: trace
!integer :: num_threads,OMP_GET_NUM_THREADS, omp_get_thread_num, numthread
integer :: iat, jat, kat, l, nat_sphere, llat
integer :: i, iiat
! double precision :: omp_get_wtime()
double precision :: omp_t1, omp_t2, omp_t3, omp_t4, omp_t5

!Variables for multielement:
integer :: num_diff_ele, iele, at_pos_counter, nat_e_cur
character(len=2), allocatable :: ele_names(:) !Short form of element name
integer, allocatable :: nat_ele(:)  !number atoms for each element
integer, allocatable :: index_ele_nat(:, :)  ! index_ele_nat(iele,iat) gives back index of original atom list(with nat at) from shortend atom list
!                                       corresponding to atom type iele with shortend position iat (e.g. index_ele_nat(2('C'),1)=16)
!  integer, allocatable :: index_nat_ele(:,:) !index_nat_ele(iele,nat) index mapping from nat_list to ele list
real*8, allocatable ::  fpall_ele(:, :)   ! All the fingerprints for all the atoms (e.g. index_nat_ele(2,16)=1)
real*8, allocatable ::  dFPdrAll_ele(:, :, :, :) !dFPdrAll(iat,ixyz,jat,iatsphere)  derivatice of the jat-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsohere, with respect to the component ixyz of atom iat
real*8, allocatable ::  dFPdalatAll_ele(:, :, :, :) !dFPdrAll(iatsphere,ixyz,jlat,l)  derivatice of the l-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of lattice vector jlat
real*8, allocatable ::  DimMatrix_ele(:, :) ! dimensionality matrix nat_ele*nat_ele
real*8, allocatable ::  dDimMdr_ele(:, :, :, :)   !dDimMdr(iat,ixyz,lat,jat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat    dDimMdalat_ele
real*8, allocatable ::  dDimMdalat_ele(:, :, :, :)   !dDimMdr(iat,ixyz,3,jat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
real*8, allocatable ::  EigenvecDimM_e(:, :)
real*8, allocatable ::  evalsDimM_e(:)   ! eigenvalues of dimensionality matrix
integer :: devalNr_e
real*8, allocatable :: dlamdadr_e(:, :, :)! dlamdadr(ixyz,iat,lat_e) is the derivative of Eval lat from DimensionalityMatrix with respect
!                                                 to all atomic positions
double precision, allocatable :: dlamdadalat_e(:, :, :)! dlamdadalat_e(ixyz,i_lattice,lat_e) is the derivative of Eval lat from DimensionalityMatrix with respect
!                                                 to latice vectors
double precision :: penalty_e, dpenaldr_e(3, nat)
double precision :: dpenaldalat_e(3,3)

integer :: ixyzmax
integer :: lat, num_threads, numthread
real*8 :: t1, t2, t3, xl, yl, zl
real*8, external :: ddot !! external function from blas package.

integer :: LDA, LWMAX, INFO, LWORK, devalNr, dwithTrace

real*8, allocatable :: fp_transpose(:,:) !fingerprint for each environment inverse

!double precision :: m_parameter_inv !Parameter for sphere size of inv_fp_calculation
double precision :: frac_num_target_atoms !number of atoms that the algorithem trys to aim for in the sphere
real*8 :: radius_cutoff

num_dlamda = num_cat !number of i for lamda_i/dR

LDA = 2*nat
LWMAX = 10000

!Which Eval should be calculated + with or without d[trace(A)-Lamda_i]/dR
devalNr = nat !Nach welchem Eigenvalue abgeleitet wird (nat= größter Eval)
dwithTrace = 1  !If dwithTrace=1 => d[trace(A)-Lamda_i]/dR || If dwithTrace=0 => only dLamda/dR

!calculate/get needed parameters:
!Calc ixyzmax
call calc_ixyzmax(width_cutoff, alat, ixyzmax)


radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff  ! radius where the polynomial is zero

do i = 1, nat
  call sym2rcov(atomnames(i), rcov(i))
  !write(*,*)atomnames(i),rcov(i)
end do
!calculate/get needed parameters finished

! write(*,*) "nat",nat
! write(*,*) '-----------------------'
! write(*,*) "alat_sym",alat
! write(*,*) '-----------------------'
! write(*,*) "rxyz_sym",rxyz
! write(*,*) '-----------------------'
! write(*,*) "atomnames",atomnames
! write(*,*) '-----------------------'
! write(*,*) "natx_sphere, ns, np",natx_sphere, ns, np
! write(*,*) '-----------------------'
! write(*,*) "width_cutoff, num_element_enviroments",width_cutoff, num_cat
! write(*,*) '-----------------------'
! write(*,*) "nex_cutoff, ixyzmax",nex_cutoff, ixyzmax
! write(*,*) '-----------------------'
! write(*,*) "rcov",rcov
! write(*,*) '-----------------------'
! write(*,*) "lengthfp",lengthfp
! write(*,*) '-----------------------'

! !Calculate multielement components:
! k=1
! do iat=1,nat
!   if (any(atomnames(iat) == ))
! end do
call get_num_ele(nat, atomnames, num_diff_ele)
!  write(*,*) "num_diff_ele",num_diff_ele

allocate (ele_names(num_diff_ele), nat_ele(num_diff_ele))
allocate (index_ele_nat(num_diff_ele, nat))
call get_element_names(nat, num_diff_ele, atomnames, ele_names)
!write (*, *) "ele_names: ", ele_names(1), " , ", ele_names(2)

!Get num_atoms per element
do iele = 1, num_diff_ele
  nat_ele(iele) = count(ele_names(iele) == atomnames)
end do
!write (*, *) "nat_ele", nat_ele

do iele = 1, num_diff_ele
  at_pos_counter = 1
  do iat = 1, nat
   if (atomnames(iat) == ele_names(iele)) then
     index_ele_nat(iele, at_pos_counter) = iat
     at_pos_counter = at_pos_counter + 1
   end if
  end do
end do
!  WRITE(*,*)"index_ele",index_ele(2,:)

!write(*,*) "shape",shape(dDimMdr_ele)

if (natx_sphere*(ns + np*3) .gt. lengthfp) stop 'increase lengthfp'

allocate (rxyz_sphere(3, natx_sphere), rcov_sphere(natx_sphere))
allocate (amplitude(natx_sphere), deramplitude(natx_sphere))
allocate (fp_loc(natx_sphere*(ns + np*3)), dfpdr(3, natx_sphere, natx_sphere*(ns + np*3)))
allocate (dfpdr0(3, natx_sphere, natx_sphere*(ns + np*3)))
allocate (dFPdrAll(nat, 3, nat, natx_sphere*(ns + np*3)), dDimMdr(nat, 3, nat, nat)) !matrixA nat times nat matrix with A_ij being an 3 times 32 object
allocate (fp((ns + 3*np)*natx_sphere, nat), DimMatrix(nat, nat), evalsDimM(nat))
allocate (WORK(LWMAX), EigenvecDimM(nat, nat), tmpResdL(nat, 3, nat))
allocate (indat(natx_sphere), nat_sphere_array(nat))
allocate (fp_transpose(nat, natx_sphere*(ns+3*np)))   !fingerprint for each environment inversed
allocate (dfp_marco(natx_sphere*(ns + np*3), 3, nat, nat))
allocate (dFPdalatAll(nat, 3, 3, natx_sphere*(ns + np*3)))

! REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np)), INTENT(OUT) :: fp  !fingerprint vectors
! REAL(8), DIMENSION(nat,3,nat,natx_sphere*(ns+3*np)), INTENT(OUT) :: dfp !fingerprint derivatives

dDimMdr = 0.d0
dFPdrAll = 0.d0
fp = 0.d0
dfpdr0 = 0.d0
dfpdr = 0.d0

!call cpu_time(t1)

  !Calculate "classic" overlap-fp and derivatives
  call atomic_fp(nat, alat, rxyz, atomnames, natx_sphere, ns, np, width_cutoff, nex_cutoff, &
                        ixyzmax, rcov, lengthfp, nat_sphere_current_max, fp, dFPdrAll, dFPdalatAll)
! call cpu_time(t2)
! print '("Time (For atomic_fp total)= ",f16.3," seconds.")',t2-t1
!write(*,*) "nat_sphere_current_max",nat_sphere_current_max,natx_sphere

penalty = 0.d0
dpenaldr = 0.d0
dpenaldalat = 0.d0

do iele = 1, num_diff_ele
  !Here start per Atom calculation of SymmetryPenalty:
  nat_e_cur = nat_ele(iele)

  !If only one atom for an element no symmetry calculation can be performed
  if ( nat_e_cur == 1) then
   cycle
  end if
  !allocate run spezific multi_ele variables:

  allocate (fpall_ele((ns + 3*np)*natx_sphere, nat_e_cur))
  allocate (dFPdrAll_ele(nat_e_cur, 3, nat, natx_sphere*(ns + np*3)))
  allocate (DimMatrix_ele(nat_e_cur, nat_e_cur))
  allocate (dDimMdr_ele(nat_e_cur, 3, nat, nat_e_cur))
  allocate (EigenvecDimM_e(nat_e_cur, nat_e_cur))
  allocate (evalsDimM_e(nat_e_cur))
  allocate (dlamdadr_e(3, nat, nat_e_cur))
  allocate (dFPdalatAll_ele(nat_e_cur, 3, 3, natx_sphere*(ns + np*3)))
  allocate (dDimMdalat_ele(nat_e_cur, 3, 3, nat_e_cur))
  allocate (dlamdadalat_e(3, 3, nat_e_cur))


  fpall_ele = 0.d0
  dFPdrAll_ele = 0.d0
  DimMatrix_ele = 0.d0
  dDimMdr_ele = 0.d0
  EigenvecDimM_e = 0.d0
  evalsDimM_e = 0.d0
  penalty_e = 0.d0
  dpenaldr_e = 0.d0
  dpenaldalat_e = 0.d0

  !num of highest eval (is last evalue)
  devalNr_e = nat_e_cur

  !Map fp/dFPdrAll onto fpall_ele/dFPdrAll_ele
  do iat = 1, nat_ele(iele)
   iiat = index_ele_nat(iele, iat)

   fpall_ele(:, iat) = fp(:, iiat)
   do jat = 1, nat
     do l = 1, nat_sphere_current_max*(ns + 3*np)
       dFPdrAll_ele(iat, 1, jat, l) = dFPdrAll(iiat, 1, jat, l)
       dFPdrAll_ele(iat, 2, jat, l) = dFPdrAll(iiat, 2, jat, l)
       dFPdrAll_ele(iat, 3, jat, l) = dFPdrAll(iiat, 3, jat, l)

       if (jat <= 3) then !alat only having 3x3 entries
         dFPdalatAll_ele(iat, 1, jat, l) = dFPdalatAll(iiat, 1, jat, l)
         dFPdalatAll_ele(iat, 2, jat, l) = dFPdalatAll(iiat, 2, jat, l)
         dFPdalatAll_ele(iat, 3, jat, l) = dFPdalatAll(iiat, 3, jat, l)
       end if
     end do
   end do
  end do

  !call cpu_time(t2)
  !!$  omp_t2 = omp_get_wtime()

  !write(*,*)"fp_timing:", omp_t2-omp_t1

  !Calc DimensionalityMatrix_ele

  do jat = 1, nat_e_cur
   do iat = 1, nat_e_cur
     DimMatrix_ele(iat, jat) = ddot((ns + 3*np)*natx_sphere, fpall_ele(1, iat), 1, fpall_ele(1, jat), 1)
   end do
  end do

  !!$  omp_t4 = omp_get_wtime()
  call calc_dDimMdr_e(nat, nat_e_cur, ns, np, natx_sphere, fpall_ele, dFPdrAll_ele, dDimMdr_ele)
  !!$  omp_t5 = omp_get_wtime()

  !Calc dDimM/dalat_element: nat = 3
  call calc_dDimMdr_e(3, nat_e_cur, ns, np, natx_sphere, fpall_ele, dFPdalatAll_ele, dDimMdalat_ele)

  EigenvecDimM_e = DimMatrix_ele
  LWORK = -1
  !dsyev if "N"=> Eigenvalues only if "V"&INFO=0 Eigenvalues and Eigenvectors
  call dsyev('V', 'U', nat_e_cur, EigenvecDimM_e, nat_e_cur, evalsDimM_e, WORK, LWORK, INFO)
  LWORK = MIN(LWMAX, INT(WORK(1)))
  call dsyev('V', 'U', nat_e_cur, EigenvecDimM_e, nat_e_cur, evalsDimM_e, WORK, LWORK, INFO)

  trace = 0.d0
  do i = 1, nat_e_cur
   trace = trace + DimMatrix_ele(i, i)
  end do

  !Penalty=Trace-Lamda(1)-Lamda(2) ...
  penalty_e = trace
  !print*, "devalNr_e: ", devalNr_e, "num_cat -1: ", num_cat - 1
  do i = 0, num_cat - 1
   penalty_e = penalty_e - evalsDimM_e(devalNr_e - i)
  end do


  dlamdadr_e = 0.d0
  call calc_multi_dEvalDimMdr_e(nat, nat_e_cur, num_dlamda, EigenvecDimM_e, dDimMdr_ele, dlamdadr_e)

  !Calculate derivate of Evalue with respect to alat:
  dlamdadalat_e = 0.d0
  call calc_multi_dEvalDimMdr_e(3, nat_e_cur, num_dlamda, EigenvecDimM_e, dDimMdalat_ele, dlamdadalat_e)
  ! !In the moment only one Eval derivate of DimMatrix is needed for penalty
  ! dpenaldr=dlamdadr(:,:,devalNr)

  do i = 0, num_cat - 1
   do iat = 1, nat
     dpenaldr_e(1, iat) = dpenaldr_e(1, iat) + dlamdadr_e(1, iat, devalNr_e - i)
     dpenaldr_e(2, iat) = dpenaldr_e(2, iat) + dlamdadr_e(2, iat, devalNr_e - i)
     dpenaldr_e(3, iat) = dpenaldr_e(3, iat) + dlamdadr_e(3, iat, devalNr_e - i)

     if (iat <= 3) then !Derivative with respect to lattice vectors (only 3 vectors)
       dpenaldalat_e(1, iat) = dpenaldalat_e(1, iat) + dlamdadalat_e(1, iat, devalNr_e - i)
       dpenaldalat_e(2, iat) = dpenaldalat_e(2, iat) + dlamdadalat_e(2, iat, devalNr_e - i)
       dpenaldalat_e(3, iat) = dpenaldalat_e(3, iat) + dlamdadalat_e(3, iat, devalNr_e - i)
     end if
   end do
  end do

  !For calculation of d[trace(A)-Lamda_i]/dR = dA_11/dR + dA_22/dR +... -dLamda/dR
  if (dwithTrace == 1) then
     dpenaldr_e = dpenaldr_e * (-1.d0)
     dpenaldalat_e = dpenaldalat_e * (-1.d0)
     do iat = 1, nat
       do jat = 1, nat_e_cur
         dpenaldr_e(1, iat) = dpenaldr_e(1, iat) + dDimMdr_ele(jat, 1, iat, jat)
         dpenaldr_e(2, iat) = dpenaldr_e(2, iat) + dDimMdr_ele(jat, 2, iat, jat)
         dpenaldr_e(3, iat) = dpenaldr_e(3, iat) + dDimMdr_ele(jat, 3, iat, jat)

         if (iat <= 3) then
           dpenaldalat_e(1, iat) = dpenaldalat_e(1, iat) + dDimMdalat_ele(jat, 1, iat, jat)
           dpenaldalat_e(2, iat) = dpenaldalat_e(2, iat) + dDimMdalat_ele(jat, 2, iat, jat)
           dpenaldalat_e(3, iat) = dpenaldalat_e(3, iat) + dDimMdalat_ele(jat, 3, iat, jat)
         end if
       end do
     end do
  end if




  !call cpu_time(t3)
  !!$  omp_t3 = omp_get_wtime()

  ! write(*,'(3(A,e10.3))')"OMP: T1: ",omp_t2-omp_t1,"      T2: ",omp_t3-omp_t2, "      Tall ",omp_t3-omp_t1, &
  ! "      TdDimdr ",omp_t5-omp_t4

  !write(*,'(3(A,e10.3))')"t2-t1",t2-t1,"      t3-t2 ",t3-t2, "      t_all ",t3-t1

  !Add penalty_e/dpendr_e to total penalty/dpendr_e
  penalty = penalty + penalty_e
  dpenaldr = dpenaldr + dpenaldr_e
  dpenaldalat = dpenaldalat + dpenaldalat_e

  !write(*,*)"pen", penalty,penalty_e

  !deallocate element spezific variables
  deallocate (fpall_ele, dFPdrAll_ele, DimMatrix_ele, dDimMdr_ele)
  deallocate (EigenvecDimM_e, evalsDimM_e, dlamdadr_e)
  deallocate (dDimMdalat_ele, dFPdalatAll_ele, dlamdadalat_e)

end do !End do over different elements for symmetry penalty


!deallocate rest:
deallocate (tmpResdL, EigenvecDimM, dDimMdr)
deallocate (DimMatrix, fp, dFPdrAll)

deallocate (rxyz_sphere, rcov_sphere)
deallocate (amplitude, deramplitude)
deallocate (fp_loc, dfpdr, evalsDimM)
deallocate (dfpdr0)
deallocate (indat)


end subroutine SymmetryPenalty_v3

subroutine atomic_fp(nat, alat, rxyz, atomnames, natx_sphere, ns, np, width_cutoff, nex_cutoff, &
                      ixyzmax, rcov, lengthfp, nat_sphere_current_max, fp_all, dFPdrAll, dFPdalatAll)
  implicit none
  integer, intent(in) :: natx_sphere! = 35  ! maximum number of atoms in the sphere
  integer, intent(in) :: ns, np !  = 1 = 1
  integer, intent(in) :: nat
  real*8, dimension(3, nat) :: rxyz !atomic input positions (unchange on output)
  real*8, dimension(3, nat) :: rxyz_2,rxyz_3 !atomic input positions (unchange on output)
  real*8 :: xyzred(3,nat) !reduced atomic coordinates
  real*8, dimension(3, 3) :: alat ! three column vectors specifying the cell vectors
  real*8, dimension(nat) :: rcov ! covalent radii neede for fingerprints
  real*8, intent(out) ::  fp_all((ns + 3*np)*natx_sphere, nat)   ! All the fingerprints for all the atoms
  real*8, intent(out) ::  dFPdrAll(nat, 3, nat, natx_sphere*(ns + np*3)) !dFPdrAll(iatsphere,ixyz,jat,l)  derivatice of the l-th eigenvalue of the
  !                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of atom jat
  real*8, intent(out) ::  dFPdalatAll(nat, 3, 3, natx_sphere*(ns + np*3)) !dFPdrAll(iatsphere,ixyz,jlat,l)  derivatice of the l-th eigenvalue of the
  !                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of lattice vector jlat
  double precision :: rxyz_sphere_reduced(3, natx_sphere),rxyz_sphere_reduced_2(3, natx_sphere)
  double precision :: rxyz_sphere_reduced_3(3, natx_sphere),rxyz_sphere_reduced_4(3, natx_sphere)
  double precision :: rxyz_sphere_2(3,natx_sphere)
  real*8, allocatable ::  rxyz_sphere(:, :), rcov_sphere(:) ! positions and covalent radii of atoms in sphere
  real*8, allocatable ::  amplitude(:), deramplitude(:) ! amplitude modulationf function for Gaussians and its derivative
  real*8, allocatable ::  fp_tmp(:) ! eigenvlues of the overlap matric which give the finperprint vector for the atomic environment
  real*8, allocatable ::  dfpdr(:, :, :) ! dfpdr(ixyz,iat_sp,l)  derivative of l-th eigenvalues with respect
  !                                       to the component ixyz  of the position of atom iat_sp in the sphere
  real*8, allocatable ::  dfpdr0(:, :, :) ! Work array (same as dfpdr) but for a number of atoms in the sphere that is less or
  !                                           equal to natx_sphere
  integer, allocatable ::  indat(:) !contains  the indices of the atoms in the sphere
  integer, allocatable ::  nat_sphere_array(:) ! needed to calculate maximum number of atoms in sphere
  integer :: nat_sphere_current_max
  integer :: iat, jat, kat, l, nex_cutoff, nat_sphere, llat, ixyz
  integer :: i, iiat
  !! array containg the chemical symbols of each atom
  character(len=2) :: atomnames(nat)
  real*8 :: width_cutoff
  integer :: ixyzmax
  integer :: lengthfp
  integer :: lat, num_threads, numthread
  real*8 :: t1, t2, t3, xl, yl, zl

  !integer :: LDA, LWMAX, INFO, LWORK, devalNr, dwithTrace



  if (natx_sphere*(ns + np*3) .gt. lengthfp) stop 'increase lengthfp'

  allocate (rxyz_sphere(3, natx_sphere), rcov_sphere(natx_sphere))
  allocate (amplitude(natx_sphere), deramplitude(natx_sphere))
  allocate (fp_tmp(natx_sphere*(ns + np*3)), dfpdr(3, natx_sphere, natx_sphere*(ns + np*3)))
  allocate (dfpdr0(3, natx_sphere, natx_sphere*(ns + np*3)))
  allocate (indat(natx_sphere), nat_sphere_array(nat))

  ! REAL(8), DIMENSION(nat, nat_sphere_max*(ns+3*np)), INTENT(OUT) :: fp  !fingerprint vectors
  ! REAL(8), DIMENSION(nat,3,nat,natx_sphere*(ns+3*np)), INTENT(OUT) :: dfp !fingerprint derivatives



  fp_all = 0.d0
  dFPdrAll = 0.d0
  dFPdalatAll = 0.d0
  dfpdr0 = 0.d0
  dfpdr = 0.d0
  rxyz_sphere_reduced = 0.d0


  !Calculate atomic fp and derivatives:

  !!$ omp_t1 = omp_get_wtime()
  !Parallize
  numthread = 0
  num_threads = 0
  !$omp parallel private(lat,llat,indat,rxyz_sphere,&
  !$omp                  amplitude, &
  !$omp                  deramplitude,rcov_sphere,&
  !$omp                  nat_sphere,i,fp_tmp,xl,yl,zl,&
  !$omp                  dfpdr0,dfpdr,kat,iiat,&
  !$omp                  rxyz_sphere_reduced,jat,&
  !$omp                  l,numthread)

  ! num_threads=OMP_GET_NUM_THREADS()
  ! numthread = omp_get_thread_num()

  !write(*,*) num_threads, numthread !if (numthread==0)

  !$omp do schedule(static)

     !Iterate over all atoms in Cell----------------------------------------------------------------------------------------------------------
     do lat = 1, nat

     ! set fp_tmp to zero
     do i = 1, natx_sphere*(ns + 3*np)
      fp_tmp(i) = 0.d0
     end do
     dfpdr = 0.d0
     dfpdr0 = 0.d0
     rxyz_sphere_reduced=0.d0
     rxyz_sphere = 0.d0
     nat_sphere = 0
     !Calc atoms in sphere
     call atoms_sphere(width_cutoff, nex_cutoff, lat, llat, ixyzmax, nat, natx_sphere, nat_sphere, alat, rxyz, rxyz_sphere, &
                       rcov, rcov_sphere, indat, amplitude, deramplitude)
     !call atoms_sphere_v2(width_cutoff, nex_cutoff, lat, llat, ixyzmax, nat, natx_sphere, nat_sphere, alat, rxyz, rxyz_sphere, &
     !                   rxyz_sphere_reduced, rcov, rcov_sphere, indat, amplitude, deramplitude)
     !rxyz_sphere_2 = rxyz_sphere
     call cart2frac(natx_sphere, alat, rxyz_sphere, rxyz_sphere_reduced)
     !Check with multiplaying by alat
     !call frac2cart(nat, alat, rxyz_sphere_reduced, rxyz_2)
     !call frac2cart(nat, alat, rxyz_sphere_reduced_2, rxyz_3)

     ! if (lat == 53) then
     !   do i=1, nat_sphere
     !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(1,i),rxyz_sphere_reduced_2(1,i),rxyz_sphere(1,i),rxyz_sphere_2(1,i),rxyz_2(1,i),rxyz_3(1,i),"1",lat,i, "xyzred"
     !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(2,i),rxyz_sphere_reduced_2(2,i),rxyz_sphere(2,i),rxyz_sphere_2(2,i),rxyz_2(2,i),rxyz_3(2,i),"2",lat,i, "xyzred"
     !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(3,i),rxyz_sphere_reduced_2(3,i),rxyz_sphere(3,i),rxyz_sphere_2(3,i),rxyz_2(3,i),rxyz_3(3,i),"3",lat,i, "xyzred"
     !   end do
     ! end if
    ! if (lat == 53) then
    !   do i=1, nat_sphere
    !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(1,i),rxyz_sphere(1,i),rxyz_2(1,i),"1",lat,i, "xyzred"
    !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(2,i),rxyz_sphere(2,i),rxyz_2(1,i),"2",lat,i, "xyzred"
    !     write(*,'(6(es10.3,1x)a1,1x,2(i2,1x),a)') rxyz_sphere_reduced(3,i),rxyz_sphere(3,i),rxyz_2(1,i),"3",lat,i, "xyzred"
    !   end do
    ! end if

     xl = rxyz(1, lat)
     yl = rxyz(2, lat)
     zl = rxyz(3, lat)
     !Calculate FP + derivative FP with respect to R
     !call cpu_time(t1)
     call xyz2devaldr(nat_sphere, rxyz_sphere, rcov_sphere, amplitude, deramplitude, llat, xl, yl, zl, ns, np, fp_tmp, dfpdr0)
     !call cpu_time(t2)
     !print '("Time (For xyz2devaldr total)= ",f16.7," seconds.")',t2-t1
     !stop 'stoppp123'
     call reformat_devaldr(natx_sphere, nat_sphere, ns, np, dfpdr0, dfpdr)

     !Save fp_tmp into all fp_array (fp)
     fp_all(:, lat) = fp_tmp

     !Calculate derivate FP with respect to lattice vectors
     !dFPdalatAll(nat, 3, 3, natx_sphere*(ns + np*3))
      do jat = 1, nat_sphere
        do l = 1, nat_sphere*(ns+3*np)
          dFPdalatAll(lat,1,1,l) = dFPdalatAll(lat,1,1,l) - rxyz_sphere_reduced(1,jat) * dfpdr(1, jat, l)
          dFPdalatAll(lat,1,2,l) = dFPdalatAll(lat,1,2,l) - rxyz_sphere_reduced(2,jat) * dfpdr(1, jat, l)
          dFPdalatAll(lat,1,3,l) = dFPdalatAll(lat,1,3,l) - rxyz_sphere_reduced(3,jat) * dfpdr(1, jat, l)

          dFPdalatAll(lat,2,1,l) = dFPdalatAll(lat,2,1,l) - rxyz_sphere_reduced(1,jat) * dfpdr(2, jat, l)
          dFPdalatAll(lat,2,2,l) = dFPdalatAll(lat,2,2,l) - rxyz_sphere_reduced(2,jat) * dfpdr(2, jat, l)
          dFPdalatAll(lat,2,3,l) = dFPdalatAll(lat,2,3,l) - rxyz_sphere_reduced(3,jat) * dfpdr(2, jat, l)

          dFPdalatAll(lat,3,1,l) = dFPdalatAll(lat,3,1,l) - rxyz_sphere_reduced(1,jat) * dfpdr(3, jat, l)
          dFPdalatAll(lat,3,2,l) = dFPdalatAll(lat,3,2,l) - rxyz_sphere_reduced(2,jat) * dfpdr(3, jat, l)
          dFPdalatAll(lat,3,3,l) = dFPdalatAll(lat,3,3,l) - rxyz_sphere_reduced(3,jat) * dfpdr(3, jat, l)
        end do
      end do

     !Remap atoms from sphere to Cell
     !dFPdrAll(Derivate from FP of atom lat(of ref Cell, direction(x,y,z), derivate due atom iiat, fp_tmp l) (dF_i/dr)
     do kat = 1, nat_sphere
      iiat = indat(kat)
      do l = 1, nat_sphere*(ns + 3*np)
        dFPdrAll(lat, 1, iiat, l) = dFPdrAll(lat, 1, iiat, l) + dfpdr(1, kat, l)
        dFPdrAll(lat, 2, iiat, l) = dFPdrAll(lat, 2, iiat, l) + dfpdr(2, kat, l)
        dFPdrAll(lat, 3, iiat, l) = dFPdrAll(lat, 3, iiat, l) + dfpdr(3, kat, l)
      end do
     end do

     !Save nat_sphere of each iteration in array to determine nat_sphere_max
     nat_sphere_array(lat) = nat_sphere

     !write(100,*) lat, nat_sphere,atomnames(lat)
     !write(*,*) "fp_tmp:", lat, ( fp(j, lat), j = 1, 5)

     end do !End do über Atome in Cell---------------------------------------------------------------------------
     !$omp end do
     !$omp end parallel

   !safe local fp to global variable
   !fp_symm=fp

   nat_sphere_current_max = 0
   do lat = 1, nat
     if (nat_sphere_array(lat) > nat_sphere_current_max) nat_sphere_current_max = nat_sphere_array(lat)
   end do

   !stop 'overlap_fp'

   !Calculate fingerprint derivative with respect to lattice vectors:

   ! dFPdalatAll(nat, 3, 3, natx_sphere*(ns + np*3)) !dFPdrAll(iatsphere,ixyz,jlat,l)  derivatice of the l-th eigenvalue of the
   ! !                                             overlaps mareix centered at atom iatsphere, with respect to the component ixyz of lattice vector jlat

   ! return_value = fp(3,3) !sum(fp)
   !
   ! do iat = 1, nat
   !   do ixyz = 1, 3
   !     dvaluedR(ixyz, iat) = dFPdrAll(3, ixyz, iat, 3)
   !   end do
   ! end do

  !  write(*,*) "------------------------"

end subroutine atomic_fp

!! converts cartesian coordinates rxyz to reduced coordinates xyzred
subroutine cart2frac(nat, alat, rxyz, xyzred)
  implicit none
  !! Number of Atoms
  integer, intent(in) :: nat
  !! Lattice Vectors.
  real*8, intent(in), dimension(3, 3) :: alat
  !! Position of the Atoms in cartesian coorinates.
  real*8, intent(in), dimension(3, nat) :: rxyz
  !! Position of the Atoms in reduced coordinates.
  real*8, intent(out), dimension(3, nat) :: xyzred

  !private variables
  integer :: iat
  real*8, dimension(3, 3) :: alatinv
  real*8 :: div

  div = alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
        alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
        alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1)
  div = 1.d0/div
  alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))*div
  alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))*div
  alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))*div
  alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))*div
  alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))*div
  alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))*div
  alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))*div
  alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))*div
  alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))*div

  do iat = 1, nat
    xyzred(1, iat) = alatinv(1, 1)*rxyz(1, iat) + alatinv(1, 2)*rxyz(2, iat) + alatinv(1, 3)*rxyz(3, iat)
    xyzred(2, iat) = alatinv(2, 1)*rxyz(1, iat) + alatinv(2, 2)*rxyz(2, iat) + alatinv(2, 3)*rxyz(3, iat)
    xyzred(3, iat) = alatinv(3, 1)*rxyz(1, iat) + alatinv(3, 2)*rxyz(2, iat) + alatinv(3, 3)*rxyz(3, iat)
  end do
end subroutine cart2frac

!! Converts reduced coordinates xyzred to cartesian coordinates rxyz
subroutine frac2cart(nat, alat, xyzred, rxyz)
  implicit none
  !! Number of atoms.
  integer, intent(in) :: nat
  !! Lattice Vecors
  real*8, intent(in), dimension(3, 3) :: alat
  !! Position of the atoms in reduced coordinates.
  real*8, dimension(3, nat), intent(in) :: xyzred
  !! Position of the atoms in cartesian coordinates.
  real*8, dimension(3, nat), intent(out) :: rxyz
  !private variables
  integer :: iat, i, j
  real*8 :: t
  do iat = 1, nat
    do i = 1, 3
      t = 0.d0
      do j = 1, 3
        t = t + xyzred(j, iat)*alat(i, j)
      end do
      rxyz(i, iat) = t
    end do
  end do
end subroutine frac2cart

subroutine calc_ixyzmax(width_cutoff, alat, ixyzmax)
  implicit none
  double precision, intent(in) :: width_cutoff
  double precision, dimension(3,3), intent(in) :: alat
  integer, intent(out) :: ixyzmax
  integer, parameter :: nwork=9
  double precision :: vol, alatalat(3,3), eigalat(3)
   !! Lattice vectors of the periodic cell.
  integer :: i, j
  integer :: info
  double precision :: workalat(9)

  !Calc ixyzmax
  vol=(alat(1,1)*alat(2,2)*alat(3,3)-alat(1,1)*alat(2,3)*alat(3,2)- &
    alat(1,2)*alat(2,1)*alat(3,3)+alat(1,2)*alat(2,3)*alat(3,1)+ &
    alat(1,3)*alat(2,1)*alat(3,2)-alat(1,3)*alat(2,2)*alat(3,1))
  if (vol.eq.0.d0 ) then ! no periodic boundary condition
    ixyzmax=0
  else  ! periodic boundary conditions
    do i=1,3
      do j=1,3
        alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
      enddo
    enddo
    call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
    !   write(*,*) !  'alat !  EVals',eigalat !
    !write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
    ! ixyzmax determines over how many periodiv images one has to
    !search to fill the sphere with atoms
    ixyzmax= int(sqrt(1.d0/eigalat(1))*width_cutoff) + 1
    endif

end subroutine calc_ixyzmax


subroutine element_number_2_atomname(atom_number, atomname)
  implicit none
  integer, intent(in) :: atom_number
  character(len=2), intent(out) :: atomname


  select case (atom_number)
    ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
    case (1)
      atomname = 'H'
    case (2)
      atomname = 'He'
    case (3)
      atomname = 'Li'
    case (4)
      atomname = 'Be'
    case (5)
      atomname = 'B'
    case (6)
      atomname = 'C'
    case (7)
      atomname = 'N'
    case (8)
      atomname = 'O'
    case (9)
      atomname = 'F'
    case (10)
      atomname = 'Ne'
    case (11)
      atomname = 'Na'
    case (12)
      atomname = 'Mg'
    case (13)
      atomname = 'Al'
    case (14)
      atomname = 'Si'
    case (15)
      atomname = 'P'
    case (16)
      atomname = 'S'
    case (17)
      atomname = 'Cl'
    case (18)
      atomname = 'Ar'
    case (19)
      atomname = 'K'
    case (20)
      atomname = 'Ca'
    case (21)
      atomname = 'Sc'
    case (22)
      atomname = 'Ti'
    case (23)
      atomname = 'V'
    case (24)
      atomname = 'Cr'
    case (25)
      atomname = 'Mn'
    case (26)
      atomname = 'Fe'
    case (27)
      atomname = 'Co'
    case (28)
      atomname = 'Ni'
    case (29)
      atomname = 'Cu'
    case (30)
      atomname = 'Zn'
    case (31)
      atomname = 'Ga'
    case (32)
      atomname = 'Ge'
    case (33)
      atomname = 'As'
    case (34)
      atomname = 'Se'
    case (35)
      atomname = 'Br'
    case (36)
      atomname = 'Kr'
    case (37)
      atomname = 'Rb'
    case (38)
      atomname = 'Sr'
    case (39)
      atomname = 'Y'
    case (40)
      atomname = 'Zr'
    case (41)
      atomname = 'Nb'
    case (42)
      atomname = 'Mo'
    case (43)
      atomname = 'Tc'
    case (44)
      atomname = 'Ru'
    case (45)
      atomname = 'Rh'
    case (46)
      atomname = 'Pd'
    case (47)
      atomname = 'Ag'
    case (48)
      atomname = 'Cd'
    case (49)
      atomname = 'In'
    case (50)
      atomname = 'Sn'
    case (51)
      atomname = 'Sb'
    case (52)
      atomname = 'Te'
    case (53)
      atomname = 'I'
    case (54)
      atomname = 'Xe'
    ! case ('Cs')
    !   rcov = 2.25d0
    ! case ('Ba')
    !   rcov = 1.98d0
    ! case ('La')
    !   rcov = 1.69d0
      !     case('Ce')
      !     case('Pr')
      !     case('Nd')
      !     case('Pm')
      !     case('Sm')
      !     case('Eu')
      !     case('Gd')
      !     case('Tb')
      !     case('Dy')
      !     case('Ho')
      !     case('Er')
      !     case('Tm')
      !     case('Yb')
    ! case ('Lu')
    !   rcov = 1.60d0
    ! case ('Hf')
    !   rcov = 1.50d0
    ! case ('Ta')
    !   rcov = 1.38d0
    ! case ('W')
    !   rcov = 1.46d0
    ! case ('Re')
    !   rcov = 1.59d0
    ! case ('Os')
    !   rcov = 1.28d0
    ! case ('Ir')
    !   rcov = 1.37d0
    ! case ('Pt')
    !   rcov = 1.28d0
    ! case ('Au')
    !   rcov = 1.44d0
    ! case ('Hg')
    !   rcov = 1.49d0
    ! case ('Tl')
    !   rcov = 1.48d0
    ! case ('Pb')
    !   rcov = 1.47d0
    ! case ('Bi')
    !   rcov = 1.46d0
      !     case('Po')
      !     case('At')
    ! case ('Rn')
    !   rcov = 1.45d0
    ! case ('LJ')   ! Lennard Jones atom
    !   rcov = 0.25d0   ! Assuming distances are about 1
    ! case ('LA')   ! Lennard Jones atom
    !   rcov = 1.122462048309373d0
    ! case ('LB')  ! Lennard Jones atom
    !   rcov = 0.9877666025122482d0
      !     case('Fr')
      !     case('Ra')
      !     case('Ac')
      !     case('Th')
      !     case('Pa')
      !     case('U')
      !     case('Np')
      !     case('Pu')
      !     case('Am')
      !     case('Cm')
    case default
      print *, " Not recognized element number ", atom_number; stop
  end select

end subroutine


subroutine sym2rcov(sym, rcov)
! returns the covalent radius of atom with chemical symbol sym
  implicit none
  real*8  :: rcov
  character(len=2) :: sym  ! chemical symbol
  select case (trim(sym))
    ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
  case ('H')
    rcov = 0.37d0
  case ('He')
    rcov = 0.32d0
  case ('Li')
    rcov = 1.34d0
  case ('Be')
    rcov = 0.90d0
  case ('B')
    rcov = 0.82d0
  case ('C')
    rcov = 0.77d0
  case ('N')
    rcov = 0.75d0
  case ('O')
    rcov = 0.73d0
  case ('F')
    rcov = 0.71d0
  case ('Ne')
    rcov = 0.69d0
  case ('Na')
    rcov = 1.54d0
  case ('Mg')
    rcov = 1.30d0
  case ('Al')
    rcov = 1.18d0
  case ('Si')
    rcov = 1.11d0
  case ('P')
    rcov = 1.06d0
  case ('S')
    rcov = 1.02d0
  case ('Cl')
    rcov = 0.99d0
  case ('Ar')
    rcov = 0.97d0
  case ('K')
    rcov = 1.96d0
  case ('Ca')
    rcov = 1.74d0
  case ('Sc')
    rcov = 1.44d0
  case ('Ti')
    rcov = 1.36d0
  case ('V')
    rcov = 1.25d0
  case ('Cr')
    rcov = 1.27d0
  case ('Mn')
    rcov = 1.39d0
  case ('Fe')
    rcov = 1.25d0
  case ('Co')
    rcov = 1.26d0
  case ('Ni')
    rcov = 1.21d0
  case ('Cu')
    rcov = 1.38d0
  case ('Zn')
    rcov = 1.31d0
  case ('Ga')
    rcov = 1.26d0
  case ('Ge')
    rcov = 1.22d0
  case ('As')
    rcov = 1.19d0
  case ('Se')
    rcov = 1.16d0
  case ('Br')
    rcov = 1.14d0
  case ('Kr')
    rcov = 1.10d0
  case ('Rb')
    rcov = 2.11d0
  case ('Sr')
    rcov = 1.92d0
  case ('Y')
    rcov = 1.62d0
  case ('Zr')
    rcov = 1.48d0
  case ('Nb')
    rcov = 1.37d0
  case ('Mo')
    rcov = 1.45d0
  case ('Tc')
    rcov = 1.56d0
  case ('Ru')
    rcov = 1.26d0
  case ('Rh')
    rcov = 1.35d0
  case ('Pd')
    rcov = 1.31d0
  case ('Ag')
    rcov = 1.53d0
  case ('Cd')
    rcov = 1.48d0
  case ('In')
    rcov = 1.44d0
  case ('Sn')
    rcov = 1.41d0
  case ('Sb')
    rcov = 1.38d0
  case ('Te')
    rcov = 1.35d0
  case ('I')
    rcov = 1.33d0
  case ('Xe')
    rcov = 1.30d0
  case ('Cs')
    rcov = 2.25d0
  case ('Ba')
    rcov = 1.98d0
  case ('La')
    rcov = 1.69d0
    !     case('Ce')
    !     case('Pr')
    !     case('Nd')
    !     case('Pm')
    !     case('Sm')
    !     case('Eu')
    !     case('Gd')
    !     case('Tb')
    !     case('Dy')
    !     case('Ho')
    !     case('Er')
    !     case('Tm')
    !     case('Yb')
  case ('Lu')
    rcov = 1.60d0
  case ('Hf')
    rcov = 1.50d0
  case ('Ta')
    rcov = 1.38d0
  case ('W')
    rcov = 1.46d0
  case ('Re')
    rcov = 1.59d0
  case ('Os')
    rcov = 1.28d0
  case ('Ir')
    rcov = 1.37d0
  case ('Pt')
    rcov = 1.28d0
  case ('Au')
    rcov = 1.44d0
  case ('Hg')
    rcov = 1.49d0
  case ('Tl')
    rcov = 1.48d0
  case ('Pb')
    rcov = 1.47d0
  case ('Bi')
    rcov = 1.46d0
    !     case('Po')
    !     case('At')
  case ('Rn')
    rcov = 1.45d0
  case ('LJ')   ! Lennard Jones atom
    rcov = 0.25d0   ! Assuming distances are about 1
  case ('LA')   ! Lennard Jones atom
    rcov = 1.122462048309373d0
  case ('LB')  ! Lennard Jones atom
    rcov = 0.9877666025122482d0
    !     case('Fr')
    !     case('Ra')
    !     case('Ac')
    !     case('Th')
    !     case('Pa')
    !     case('U')
    !     case('Np')
    !     case('Pu')
    !     case('Am')
    !     case('Cm')
  case default
    print *, " Not recognized atomic type "//sym; stop
  end select

  rcov = rcov/0.52917720859d0   ! convert to atomic units

!  write(*,*) rcov

end subroutine sym2rcov

subroutine get_num_ele(nat, atomnames, num_diff_ele)
  implicit none
  integer :: nat

  !! array containg the chemical symbols of each atom
  character(len=2) :: atomnames(nat)
  integer :: i, num, num_diff_ele
  logical, dimension(nat) :: mask

  mask = .false.

  do i = 1, size(atomnames)
    num = count(atomnames(i) == atomnames)
    if (num == 1) then
      mask(i) = .true.
    else
      !atomnames(i)==atomnames gives array there this is true
      !mask gives array there mask is true
      !if atome name is not in "mask" then all false => mask(i)=.true
      if (.not. any(atomnames(i) == atomnames .and. mask)) mask(i) = .true.
    end if
  end do
  num_diff_ele = count(mask)
end subroutine get_num_ele

subroutine get_element_names(nat, num_diff_ele, atomnames, ele_names)
  implicit none
  integer :: nat, num_diff_ele
  !! array containg the chemical symbols of each atom
  character(len=2) :: atomnames(nat)
  character(len=2), dimension(num_diff_ele) :: ele_names

  integer :: i, num
  logical, dimension(nat) :: mask

  mask = .false.

  do i = 1, size(atomnames)
    num = count(atomnames(i) == atomnames)
    if (num == 1) then
      mask(i) = .true.
    else
      if (.not. any(atomnames(i) == atomnames .and. mask)) mask(i) = .true.
    end if
  end do
  ele_names = pack(atomnames, mask)
end subroutine get_element_names

!call calc_dDimMdr_e(nat,nat_e_cur,ns,np,natx_sphere,fpall_ele,dFPdrAll_ele,dDimMdr_ele)
subroutine calc_dDimMdr_e(nat, nat_e, ns, np, natx_sphere, fp, tmpResdL, dDimMdr)
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
  implicit none
  integer :: nat, iat, jat, kat, nat_e
  integer :: ns
  integer :: np
  integer :: natx_sphere
  real*8 :: dDimMdr(nat_e, 3, nat, nat_e) !dDimMdr(iat_e,dixyz,djat,lat_e) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: dDimMdrWork(nat_e, 3, nat, nat_e) !WorkArray
  real*8 :: tmpResdL(nat_e*nat*3, (ns + 3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
  real*8 :: fp((ns + 3*np)*natx_sphere, nat_e)

  call calc_dDimMdr_MatMul_e(nat, nat_e, ns, np, natx_sphere, fp, tmpResdL, dDimMdr)
  dDimMdrWork = dDimMdr
  do kat = 1, nat
    do jat = 1, nat_e
      do iat = 1, nat_e
        dDimMdr(iat, 1, kat, jat) = dDimMdr(iat, 1, kat, jat) + dDimMdrWork(jat, 1, kat, iat)
        dDimMdr(iat, 2, kat, jat) = dDimMdr(iat, 2, kat, jat) + dDimMdrWork(jat, 2, kat, iat)
        dDimMdr(iat, 3, kat, jat) = dDimMdr(iat, 3, kat, jat) + dDimMdrWork(jat, 3, kat, iat)
      end do
    end do
  end do
end subroutine calc_dDimMdr_e

subroutine calc_dDimMdr_MatMul_e(nat, nat_e, ns, np, natx_sphere, fp, tmpResdL, Res)!dFPdrAll,dDimMdr -> tmpResdL,Res
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
  implicit none
  integer :: nat, nat_e
  integer :: ns, np, natx_sphere
  real*8 :: tmpResdL(nat_e*nat*3, (ns + 3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
  real*8 :: fp((ns + 3*np)*natx_sphere, nat_e)
  real*8 :: Res(nat*nat*3, nat)

  double precision :: alpha, beta
  INTEGER M, K, N

!tmpResdL being an work array to remap the 4 dim matrix dFPdrAll into an 2 dim matrix
!for faster computing time while doing an matrix multiplication
!tmpResdL is structured like:

! (dFPdrAll with xyz=x,lat=1)
! (dFPdrAll with xyz=y,lat=1)
! (dFPdrAll with xyz=z,lat=1)
! (dFPdrAll with xyz=x,lat=2)
! (dFPdrAll with xyz=y,lat=2)
! (dFPdrAll with xyz=z,lat=2)
  ! .
  ! .
  ! .
  ! .
! (dFPdrAll with xyz=x,lat=nat)
! (dFPdrAll with xyz=y,lat=nat)
! (dFPdrAll with xyz=z,lat=nat)
!
!With e.g. (dFPdrAll with xyz=x,lat=1) being the derivative of the fp
!With respect to the x direction and to lat=1.
!(dFPdrAll with xyz=x,lat=1) has the dimensions (nat,natx_sphere*(ns+np*3))
!tmpResdL has therefore the Dimsions (nat*nat*3,natx_sphere*(ns+np*3))

!Parameter declaration for DGEMM MatMul
  alpha = 1.d0
  beta = 0.d0
  M = nat_e*nat*3
  K = (ns + 3*np)*natx_sphere
  N = nat_e

!DGEMM for optimized matrix multiplication
!DGEMM computes real matrix Res=alpha*tmpResdL*fp+beta*Res
  CALL DGEMM('N', 'N', M, N, K, alpha, tmpResdL, M, fp, K, beta, Res, M)

! The Res Matrix is structured like:
! (dDimMdr with xyz=x,lat=1)
! (dDimMdr with xyz=y,lat=1)
! (dDimMdr with xyz=z,lat=1)
! (dDimMdr with xyz=x,lat=2)
! (dDimMdr with xyz=y,lat=2)
! (dDimMdr with xyz=z,lat=2)
  ! .
  ! .
  ! .
  ! .
! (dDimMdr with xyz=x,lat=nat)
! (dDimMdr with xyz=y,lat=nat)
! (dDimMdr with xyz=z,lat=nat)
!
!With e.g. (dDimMdr with xyz=x,lat=1) being the derivative of the DimensionalityMatrix
!With respect to the x direction and to lat=1 (with dimesionions (nat,nat))
!Res has therefore the Dimsions (nat*nat*3,nat)
end subroutine calc_dDimMdr_MatMul_e

subroutine calc_dDimMdr(nat, ns, np, natx_sphere, fp, tmpResdL, dDimMdr)
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
  implicit none
  integer :: nat, iat, jat, kat
  integer :: ns, np, natx_sphere
  real*8 :: dDimMdr(nat, 3, nat, nat) !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: dDimMdrWork(nat, 3, nat, nat) !WorkArray
  real*8 :: tmpResdL(nat*nat*3, (ns + 3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
  real*8 :: fp((ns + 3*np)*natx_sphere, nat)

  call calc_dDimMdr_MatMul(nat, ns, np, natx_sphere, fp, tmpResdL, dDimMdr)
  dDimMdrWork = dDimMdr
  do kat = 1, nat
    do jat = 1, nat
      do iat = 1, nat
        dDimMdr(iat, 1, kat, jat) = dDimMdr(iat, 1, kat, jat) + dDimMdrWork(jat, 1, kat, iat)
        dDimMdr(iat, 2, kat, jat) = dDimMdr(iat, 2, kat, jat) + dDimMdrWork(jat, 2, kat, iat)
        dDimMdr(iat, 3, kat, jat) = dDimMdr(iat, 3, kat, jat) + dDimMdrWork(jat, 3, kat, iat)
      end do
    end do
  end do
end subroutine calc_dDimMdr

subroutine calc_dDimMdr_MatMul(nat, ns, np, natx_sphere, fp, tmpResdL, Res)!dFPdrAll,dDimMdr -> tmpResdL,Res
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
  implicit none
  integer :: nat
  integer :: ns, np, natx_sphere
  real*8 :: tmpResdL(nat*nat*3, (ns + 3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
  real*8 :: fp((ns + 3*np)*natx_sphere, nat)
  real*8 :: Res(nat*nat*3, nat)

  double precision :: alpha, beta
  INTEGER M, K, N

!tmpResdL being an work array to remap the 4 dim matrix dFPdrAll into an 2 dim matrix
!for faster computing time while doing an matrix multiplication
!tmpResdL is structured like:

! (dFPdrAll with xyz=x,lat=1)
! (dFPdrAll with xyz=y,lat=1)
! (dFPdrAll with xyz=z,lat=1)
! (dFPdrAll with xyz=x,lat=2)
! (dFPdrAll with xyz=y,lat=2)
! (dFPdrAll with xyz=z,lat=2)
  ! .
  ! .
  ! .
  ! .
! (dFPdrAll with xyz=x,lat=nat)
! (dFPdrAll with xyz=y,lat=nat)
! (dFPdrAll with xyz=z,lat=nat)
!
!With e.g. (dFPdrAll with xyz=x,lat=1) being the derivative of the fp
!With respect to the x direction and to lat=1.
!(dFPdrAll with xyz=x,lat=1) has the dimensions (nat,natx_sphere*(ns+np*3))
!tmpResdL has therefore the Dimsions (nat*nat*3,natx_sphere*(ns+np*3))

!Parameter declaration for DGEMM MatMul
  alpha = 1.d0
  beta = 0.d0
  M = nat*nat*3
  K = (ns + 3*np)*natx_sphere
  N = nat

!DGEMM for optimized matrix multiplication
!DGEMM computes real matrix Res=alpha*tmpResdL*fp+beta*Res
  CALL DGEMM('N', 'N', M, N, K, alpha, tmpResdL, M, fp, K, beta, Res, M)

! The Res Matrix is structured like:
! (dDimMdr with xyz=x,lat=1)
! (dDimMdr with xyz=y,lat=1)
! (dDimMdr with xyz=z,lat=1)
! (dDimMdr with xyz=x,lat=2)
! (dDimMdr with xyz=y,lat=2)
! (dDimMdr with xyz=z,lat=2)
  ! .
  ! .
  ! .
  ! .
! (dDimMdr with xyz=x,lat=nat)
! (dDimMdr with xyz=y,lat=nat)
! (dDimMdr with xyz=z,lat=nat)
!
!With e.g. (dDimMdr with xyz=x,lat=1) being the derivative of the DimensionalityMatrix
!With respect to the x direction and to lat=1 (with dimesionions (nat,nat))
!Res has therefore the Dimsions (nat*nat*3,nat)
end subroutine calc_dDimMdr_MatMul

subroutine calc_multi_dEvalDimMdr_e(nat, nat_e, num_dlamda, EigenvecDimM_e, dDimMdr_e, dlamdadr_e)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
  implicit none
  integer :: nat, nat_e, iat, num_dlamda, i
  real*8 :: EigenvecDimM_e(nat_e, nat_e) !EigenvecDimM(iat,lat) Evec of DimMatrix belonging to Eval lat for dEval/dr
  real*8 :: dDimMdr_e(nat_e, 3, nat, nat_e)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: dlamdadr_e(3, nat, nat_e) ! dlamdadr(ixyz,iat,lat) is the derivative of Eval lat from DimensionalityMatrix with respect
  !                                            to all atomic positions

  if (num_dlamda > nat_e) stop "num_dlamda>nat_e"

!$omp parallel private(iat, i)

!$omp do schedule(static)
  do i = 1, num_dlamda
!iat=Evalue that we calulate the derivate of
    iat = nat_e - i + 1

    ! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
    ! Returns dpenaldr
    call calc_dEvalDimMdr_e(nat, nat_e, EigenvecDimM_e(1, iat), dDimMdr_e, dlamdadr_e(1, 1, iat))

  end do !End do über Atome in Cell---------------------------------------------------------------------------
!$omp end do
!$omp end parallel

! Old single derivative call in SymmetryPenalty:
! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
! Returns dpenaldr
!call calc_dEvalDimMdr(nat,EigenvecDimM(:,devalNr),dDimMdr,dpenaldr)
end subroutine calc_multi_dEvalDimMdr_e

subroutine calc_dEvalDimMdr_e(nat, nat_e, evec_e, dDimMdr_e, devalDimMdr)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
  implicit none
  integer :: nat, nat_e, iat, jat, kat
  real*8 :: evec_e(nat_e) !Evec of DimMatrix belonging to Eval for dEval/dr
  real*8 :: dDimMdr_e(nat_e, 3, nat, nat_e)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: tmpResdL(nat_e, 3, nat) !Working Array
  real*8 :: devalDimMdr(3, nat) ! devalDimMdr(ixyz,lat) Derivate of Eval from DimMatrix with respect to ixyz,lat

!Calc derivate form Eval from DimMatrix to dr: dLamda_i/dr=<v_i|dDimM/dr|v_i>
  tmpResdL = 0.d0
  do iat = 1, nat_e
    do jat = 1, nat
      do kat = 1, nat_e
        tmpResdL(iat, 1, jat) = tmpResdL(iat, 1, jat) + evec_e(kat)*dDimMdr_e(kat, 1, jat, iat)
        tmpResdL(iat, 2, jat) = tmpResdL(iat, 2, jat) + evec_e(kat)*dDimMdr_e(kat, 2, jat, iat)
        tmpResdL(iat, 3, jat) = tmpResdL(iat, 3, jat) + evec_e(kat)*dDimMdr_e(kat, 3, jat, iat)
      end do
    end do
  end do

  devalDimMdr = 0.d0
  do iat = 1, nat
    do kat = 1, nat_e
      devalDimMdr(1, iat) = devalDimMdr(1, iat) + tmpResdL(kat, 1, iat)*evec_e(kat)
      devalDimMdr(2, iat) = devalDimMdr(2, iat) + tmpResdL(kat, 2, iat)*evec_e(kat)
      devalDimMdr(3, iat) = devalDimMdr(3, iat) + tmpResdL(kat, 3, iat)*evec_e(kat)
    end do
  end do

end subroutine calc_dEvalDimMdr_e

subroutine calc_multi_dEvalDimMdr(nat, num_dlamda, EigenvecDimM, dDimMdr, dlamdadr)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
  implicit none
  integer :: nat, iat, num_dlamda, i
  real*8 :: EigenvecDimM(nat, nat) !EigenvecDimM(iat,lat) Evec of DimMatrix belonging to Eval lat for dEval/dr
  real*8 :: dDimMdr(nat, 3, nat, nat)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: dlamdadr(3, nat, nat) ! dlamdadr(ixyz,iat,lat) is the derivative of Eval lat from DimensionalityMatrix with respect
  !                                            to all atomic positions

  if (num_dlamda > nat) stop "num_dlamda>nat"

!$omp parallel private(iat, i)

!$omp do schedule(static)
  do i = 1, num_dlamda
!iat=Evalue that we calulate the derivate of
    iat = nat - i + 1

    ! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
    ! Returns dpenaldr
    call calc_dEvalDimMdr(nat, EigenvecDimM(1, iat), dDimMdr, dlamdadr(1, 1, iat))

  end do !End do über Atome in Cell---------------------------------------------------------------------------
!$omp end do
!$omp end parallel

! Old single derivative call in SymmetryPenalty:
! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
! Returns dpenaldr
!call calc_dEvalDimMdr(nat,EigenvecDimM(:,devalNr),dDimMdr,dpenaldr)
end subroutine calc_multi_dEvalDimMdr

subroutine calc_dEvalDimMdr(nat, Evec, dDimMdr, devalDimMdr)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
  implicit none
  integer :: nat, iat, jat, kat
  real*8 :: Evec(nat) !Evec of DimMatrix belonging to Eval for dEval/dr
  real*8 :: dDimMdr(nat, 3, nat, nat)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8 :: tmpResdL(nat, 3, nat) !Working Array
  real*8 :: devalDimMdr(3, nat) ! devalDimMdr(ixyz,lat) Derivate of Eval from DimMatrix with respect to ixyz,lat

!Calc derivate form Eval from DimMatrix to dr: dLamda_i/dr=<v_i|dDimM/dr|v_i>
  tmpResdL = 0.d0
  do iat = 1, nat
    do jat = 1, nat
      do kat = 1, nat
        tmpResdL(iat, 1, jat) = tmpResdL(iat, 1, jat) + Evec(kat)*dDimMdr(kat, 1, jat, iat)
        tmpResdL(iat, 2, jat) = tmpResdL(iat, 2, jat) + Evec(kat)*dDimMdr(kat, 2, jat, iat)
        tmpResdL(iat, 3, jat) = tmpResdL(iat, 3, jat) + Evec(kat)*dDimMdr(kat, 3, jat, iat)
      end do
    end do
  end do

  devalDimMdr = 0.d0
  do iat = 1, nat
    do kat = 1, nat
      devalDimMdr(1, iat) = devalDimMdr(1, iat) + tmpResdL(kat, 1, iat)*Evec(kat)
      devalDimMdr(2, iat) = devalDimMdr(2, iat) + tmpResdL(kat, 2, iat)*Evec(kat)
      devalDimMdr(3, iat) = devalDimMdr(3, iat) + tmpResdL(kat, 3, iat)*Evec(kat)
    end do
  end do
end subroutine calc_dEvalDimMdr

subroutine reformat_devaldr(natx_sphere, nat_sphere, ns, np, devaldr0, devaldr)
  implicit none
  integer :: natx_sphere, nat_sphere, ns, np, i, j
  real*8, dimension(3, nat_sphere, nat_sphere*(ns + np*3)) :: devaldr0
  real*8, dimension(3, natx_sphere, natx_sphere*(ns + np*3)) :: devaldr

  do j = 1, nat_sphere*(ns + np*3)
    do i = 1, nat_sphere
      devaldr(1, i, j) = devaldr0(1, i, j)
      devaldr(2, i, j) = devaldr0(2, i, j)
      devaldr(3, i, j) = devaldr0(3, i, j)
    end do
  end do
end subroutine reformat_devaldr

subroutine xyz2devaldr(nat, rxyz, rcov, amplitude, deramplitude, lat, xl, yl, zl, ns, np, eval, devaldr_out)
! Calculqates the derivative of all eigenvalues of an atomic fingerprint with respect to the atomic positions
! Nat=NatSphere da wir eval von sphere davor hatten bzw. für mich sollte es
! Eval = FP
  implicit none
  integer :: nat
  real*8 :: rxyz(3, nat), rcov(nat), eval(nat*(ns + np*3)), devaldr_out(3, nat, nat*(ns + np*3))
  real*8 :: amplitude(nat), deramplitude(nat)
  real*8, allocatable ::  ovrlp(:, :), evecn(:, :), evec(:, :), work(:), devaldr(:, :, :)
  real*8 alpha(nat), cs(10), cp(10)
  integer :: norb
  integer :: lat
  real*8 :: xl
  real*8 :: yl
  real*8 :: zl
  integer :: ns
  integer :: np

  real*8 :: ai, aj, deri1, deri2, deri3, derj1, derj2, derj3, derl1, derl2, derl3
  real*8 :: derx00i, derx00j, derx00l, derx01i, derx01j, derx01l, derx02i, derx02j, derx02l
  real*8 :: derx10i, derx10j, derx10l, derx11i, derx11j, derx11l
  real*8 :: derx12i, derx12j, derx12l, derx20i, derx20j, derx20l
  real*8 :: derx21i, derx21j, derx21l, derx22i, derx22j, derx22l

  real*8 :: dery00i, dery00j, dery00l, dery01i, dery01j, dery01l
  real*8 :: dery02i, dery02j, dery02l, dery10i, dery10j, dery10l
  real*8 :: dery11i, dery11j, dery11l, dery12i, dery12j, dery12l
  real*8 :: dery20i, dery20j, dery20l, dery21i, dery21j, dery21l
  real*8 :: dery22i, dery22j, dery22l

  real*8 :: derz00i, derz00j, derz00l, derz01i, derz01j, derz01l
  real*8 :: derz02i, derz02j, derz02l, derz10i, derz10j, derz10l
  real*8 :: derz11i, derz11j, derz11l, derz12i, derz12j, derz12l
  real*8 :: derz20i, derz20j, derz20l, derz21i, derz21j, derz21l
  real*8 :: derz22i, derz22j, derz22l

  real*8 :: dipjampl, djpiampl
  integer :: i, iat, info, iorb, ip, is, jat, jorb, jp, js, kkorb, korb, lwork, nt
  real*8 :: pevec, pevec00, pevec01, pevec02, pevec10, pevec11, pevec12, pevec20
  real*8 :: pevec21, pevec22, pijampl
  real*8 :: r2, sij, t1, t2, t3, t4, t5, tt, xi, xij, xil, xj, xjl, yi, yij, yil
  real*8 :: yjl, zi, zij, zil, zj, zjl, yj

  real :: start, finish

  nt = 3*np + ns
  norb = nat*(ns + np*3)
  allocate (ovrlp(norb, norb), evecn(norb, norb), evec(norb, norb), devaldr(norb, 3, nat))
  lwork = 100*norb
  allocate (work(lwork))

  do iat = 1, nat
    alpha(iat) = .5d0/rcov(iat)**2
  end do
  ! Specify the width of the Gaussians if several Gaussians per l-channel are used
  do i = 1, 10
    cs(i) = sqrt(2.d0)**(i - 1)
    cp(i) = sqrt(2.d0)**(i - 1)
  end do

  call crtovrlp(nat, rxyz, alpha, cs, cp, ns, np, ovrlp)

  do iat = 1, nat
    do iorb = 1, norb
      devaldr(iorb, 1, iat) = 0.d0
      devaldr(iorb, 2, iat) = 0.d0
      devaldr(iorb, 3, iat) = 0.d0
    end do
  end do

  call multamp(nat, ovrlp, amplitude, norb, ns, np, evecn)
  call cpu_time(start)
  call dsyev('V', 'L', norb, evecn, norb, eval, work, lwork, info)
  if (info /= 0) stop ' ERROR in dsyev'
  call cpu_time(finish)
  !print '("Time (For in xyz2devaldr: dsyev)= ",f16.7," seconds.")',finish-start
  ! eigenvalues in decreasing order
  do i = 1, norb/2
    t1 = eval(i)
    t2 = eval(norb - i + 1)
    eval(i) = t2
    eval(norb - i + 1) = t1
  end do

!return    ! if no derivatives are needed
  call rots(norb, norb, evecn, evec)

! Now calculate derivatives
  !  <s|s>
  do jat = 1, nat
    do js = 1, ns
      jorb = (jat - 1)*nt + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          iorb = (iat - 1)*nt + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          ! derivatives
          tt = -2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

          pijampl = amplitude(iat)*amplitude(jat)
          dipjampl = deramplitude(iat)*amplitude(jat)
          djpiampl = deramplitude(jat)*amplitude(iat)

          deri1 = pijampl*(tt*ovrlp(iorb, jorb)*xij) + dipjampl*ovrlp(iorb, jorb)*xil
          derj1 = -pijampl*(tt*ovrlp(iorb, jorb)*xij) + djpiampl*ovrlp(iorb, jorb)*xjl
          derl1 = -xil*dipjampl*ovrlp(iorb, jorb) - djpiampl*ovrlp(iorb, jorb)*xjl

          deri2 = pijampl*(tt*ovrlp(iorb, jorb)*yij) + dipjampl*ovrlp(iorb, jorb)*yil
          derj2 = -pijampl*(tt*ovrlp(iorb, jorb)*yij) + djpiampl*ovrlp(iorb, jorb)*yjl
          derl2 = -yil*dipjampl*ovrlp(iorb, jorb) - djpiampl*ovrlp(iorb, jorb)*yjl

          deri3 = pijampl*(tt*ovrlp(iorb, jorb)*zij) + dipjampl*ovrlp(iorb, jorb)*zil
          derj3 = -pijampl*(tt*ovrlp(iorb, jorb)*zij) + djpiampl*ovrlp(iorb, jorb)*zjl
          derl3 = -zil*dipjampl*ovrlp(iorb, jorb) - djpiampl*ovrlp(iorb, jorb)*zjl

          do korb = 1, norb
            kkorb = norb - korb + 1  ! Also deriavtive in decreasing order of eigenvalues
            pevec = evec(korb, iorb)*evec(korb, jorb)

            devaldr(kkorb, 1, iat) = devaldr(kkorb, 1, iat) + pevec*deri1
            devaldr(kkorb, 1, jat) = devaldr(kkorb, 1, jat) + pevec*derj1
            devaldr(kkorb, 1, lat) = devaldr(kkorb, 1, lat) + pevec*derl1

            devaldr(kkorb, 2, iat) = devaldr(kkorb, 2, iat) + pevec*deri2
            devaldr(kkorb, 2, jat) = devaldr(kkorb, 2, jat) + pevec*derj2
            devaldr(kkorb, 2, lat) = devaldr(kkorb, 2, lat) + pevec*derl2

            devaldr(kkorb, 3, iat) = devaldr(kkorb, 3, iat) + pevec*deri3
            devaldr(kkorb, 3, jat) = devaldr(kkorb, 3, jat) + pevec*derj3
            devaldr(kkorb, 3, lat) = devaldr(kkorb, 3, lat) + pevec*derl3
          end do
        end do
      end do
    end do
  end do

  if (np .eq. 0) goto 1111

  !  <pi|sj>
  do jat = 1, nat ! kat, kat ! 1,nat
    do js = 1, ns

      jorb = (jat - 1)*nt + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np

          iorb = (iat - 1)*nt + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = -2.d0*sqrt(ai)*aj/t2

          ! derivatives
          t5 = -2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

          pijampl = amplitude(iat)*amplitude(jat)
          !dipjampl=deramplitude(iat)*amplitude(jat)
          !djpiampl=deramplitude(jat)*amplitude(iat)

          derx00i = pijampl*(ovrlp(iorb, jorb)*t5*xij + t3*sij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derx00j = -pijampl*(ovrlp(iorb, jorb)*t5*xij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl
          derx00l = -xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl

          dery00i = pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          dery00j = -pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl
          dery00l = -yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl

          derz00i = pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derz00j = -pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl
          derz00l = -zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl

          derx10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          derx10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*xjl
          derx10l = -xil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*xjl

          dery10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*yij + t3*sij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          dery10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*yij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*yjl
          dery10l = -yil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*yjl

          derz10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          derz10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*zjl
          derz10l = -zil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*zjl

          derx20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          derx20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*xjl
          derx20l = -xil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*xjl

          dery20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          dery20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*yjl
          dery20l = -yil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*yjl

          derz20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*zij + t3*sij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          derz20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*zij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*zjl
          derz20l = -zil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*zjl

          do korb = 1, norb
            kkorb = norb - korb + 1  ! Also derivative in decreasing order of eigenvalues
            pevec00 = evec(korb, iorb + 0)*evec(korb, jorb)
            pevec10 = evec(korb, iorb + 1)*evec(korb, jorb)
            pevec20 = evec(korb, iorb + 2)*evec(korb, jorb)

            devaldr(kkorb, 1, iat) = devaldr(kkorb, 1, iat) + pevec00*derx00i + pevec10*derx10i + pevec20*derx20i
            devaldr(kkorb, 1, jat) = devaldr(kkorb, 1, jat) + pevec00*derx00j + pevec10*derx10j + pevec20*derx20j
            devaldr(kkorb, 1, lat) = devaldr(kkorb, 1, lat) + pevec00*derx00l + pevec10*derx10l + pevec20*derx20l

            devaldr(kkorb, 2, iat) = devaldr(kkorb, 2, iat) + pevec00*dery00i + pevec10*dery10i + pevec20*dery20i
            devaldr(kkorb, 2, jat) = devaldr(kkorb, 2, jat) + pevec00*dery00j + pevec10*dery10j + pevec20*dery20j
            devaldr(kkorb, 2, lat) = devaldr(kkorb, 2, lat) + pevec00*dery00l + pevec10*dery10l + pevec20*dery20l

            devaldr(kkorb, 3, iat) = devaldr(kkorb, 3, iat) + pevec00*derz00i + pevec10*derz10i + pevec20*derz20i
            devaldr(kkorb, 3, jat) = devaldr(kkorb, 3, jat) + pevec00*derz00j + pevec10*derz10j + pevec20*derz20j
            devaldr(kkorb, 3, lat) = devaldr(kkorb, 3, lat) + pevec00*derz00l + pevec10*derz10l + pevec20*derz20l

          end do
        end do
      end do
    end do
  end do

  !  <si|pj>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*nt + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          iorb = (iat - 1)*nt + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = +2.d0*sqrt(aj)*ai/t2

          ! derivatives
          !tt= -2.d0*t1/t2 * sij
          t5 = -2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

          pijampl = amplitude(iat)*amplitude(jat)
          !dipjampl=deramplitude(iat)*amplitude(jat)
          !djpiampl=deramplitude(jat)*amplitude(iat)

          derx00i = pijampl*(ovrlp(iorb, jorb)*t5*xij + t3*sij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derx00j = -pijampl*(ovrlp(iorb, jorb)*t5*xij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl
          derx00l = -xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl

          dery00i = pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          dery00j = -pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl
          dery00l = -yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl

          derz00i = pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derz00j = -pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl
          derz00l = -zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl

          derx01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          derx01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*xjl
          derx01l = -xil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*xjl

          dery01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*yij + t3*sij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          dery01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*yij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*yjl
          dery01l = -yil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*yjl

          derz01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          derz01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*zjl
          derz01l = -zil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*zjl

          derx02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          derx02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*xjl
          derx02l = -xil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*xjl

          dery02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          dery02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*yjl
          dery02l = -yil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*yjl

          derz02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*zij + t3*sij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          derz02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*zij + t3*sij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*zjl
          derz02l = -zil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*zjl

          do korb = 1, norb
            kkorb = norb - korb + 1  ! Also derivative in decreasing order of eigenvalues
            pevec00 = evec(korb, iorb)*evec(korb, jorb + 0)
            pevec01 = evec(korb, iorb)*evec(korb, jorb + 1)
            pevec02 = evec(korb, iorb)*evec(korb, jorb + 2)

            devaldr(kkorb, 1, iat) = devaldr(kkorb, 1, iat) + pevec00*derx00i + pevec01*derx01i + pevec02*derx02i
            devaldr(kkorb, 1, jat) = devaldr(kkorb, 1, jat) + pevec00*derx00j + pevec01*derx01j + pevec02*derx02j
            devaldr(kkorb, 1, lat) = devaldr(kkorb, 1, lat) + pevec00*derx00l + pevec01*derx01l + pevec02*derx02l

            devaldr(kkorb, 2, iat) = devaldr(kkorb, 2, iat) + pevec00*dery00i + pevec01*dery01i + pevec02*dery02i
            devaldr(kkorb, 2, jat) = devaldr(kkorb, 2, jat) + pevec00*dery00j + pevec01*dery01j + pevec02*dery02j
            devaldr(kkorb, 2, lat) = devaldr(kkorb, 2, lat) + pevec00*dery00l + pevec01*dery01l + pevec02*dery02l

            devaldr(kkorb, 3, iat) = devaldr(kkorb, 3, iat) + pevec00*derz00i + pevec01*derz01i + pevec02*derz02i
            devaldr(kkorb, 3, jat) = devaldr(kkorb, 3, jat) + pevec00*derz00j + pevec01*derz01j + pevec02*derz02j
            devaldr(kkorb, 3, lat) = devaldr(kkorb, 3, lat) + pevec00*derz00l + pevec01*derz01l + pevec02*derz02l
          end do

        end do
      end do
    end do
  end do

  !  <p|p>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*nt + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np
          iorb = (iat - 1)*nt + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t4 = 2.d0*sqrt(t1)/t2
          t5 = -2.d0*t1/t2

          ! derivatives

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl
          pijampl = amplitude(iat)*amplitude(jat)
          !dipjampl=deramplitude(iat)*amplitude(jat)
          !djpiampl=deramplitude(jat)*amplitude(iat)

          derx00i = pijampl*(ovrlp(iorb, jorb)*t5*xij + t4*t5*sij*xij*2.d0) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derx00j = -pijampl*(ovrlp(iorb, jorb)*t5*xij + t4*t5*sij*xij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl
          derx00l = -xil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*xjl

          dery00i = pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          dery00j = -pijampl*(ovrlp(iorb, jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl
          dery00l = -yil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*yjl

          derz00i = pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat)
          derz00j = -pijampl*(ovrlp(iorb, jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl
          derz00l = -zil*deramplitude(iat)*ovrlp(iorb, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb)*deramplitude(jat)*zjl

          derx10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*xij + t4*t5*sij*yij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          derx10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*xij + t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*xjl
          derx10l = -xil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*xjl

          dery10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*yij + t4*t5*sij*xij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          dery10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*yij + t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*yjl
          dery10l = -yil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*yjl

          derz10i = pijampl*(ovrlp(iorb + 1, jorb)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat)
          derz10j = -pijampl*(ovrlp(iorb + 1, jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*zjl
          derz10l = -zil*deramplitude(iat)*ovrlp(iorb + 1, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb)*deramplitude(jat)*zjl

          derx20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*xij + t4*t5*sij*zij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          derx20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*xij + t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*xjl
          derx20l = -xil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*xjl

          dery20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          dery20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*yjl
          dery20l = -yil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*yjl

          derz20i = pijampl*(ovrlp(iorb + 2, jorb)*t5*zij + t4*t5*sij*xij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat)
          derz20j = -pijampl*(ovrlp(iorb + 2, jorb)*t5*zij + t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*zjl
          derz20l = -zil*deramplitude(iat)*ovrlp(iorb + 2, jorb)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb)*deramplitude(jat)*zjl

          derx01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*xij + t4*t5*sij*yij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          derx01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*xij + t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*xjl
          derx01l = -xil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*xjl

          dery01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*yij + t4*t5*sij*xij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          dery01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*yij + t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*yjl
          dery01l = -yil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*yjl

          derz01i = pijampl*(ovrlp(iorb, jorb + 1)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat)
          derz01j = -pijampl*(ovrlp(iorb, jorb + 1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*zjl
          derz01l = -zil*deramplitude(iat)*ovrlp(iorb, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 1)*deramplitude(jat)*zjl

          derx11i = pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat)
          derx11j = -pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*xjl
          derx11l = -xil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*xjl

          dery11i = pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*yij + t4*t5*sij*yij*2.d0) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat)
          dery11j = -pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*yij + t4*t5*sij*yij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*yjl
          dery11l = -yil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*yjl

          derz11i = pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*zij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat)
          derz11j = -pijampl*(ovrlp(iorb + 1, jorb + 1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*zjl
          derz11l = -zil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 1)*deramplitude(jat)*zjl

          derx21i = pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat)
          derx21j = -pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*xjl
          derx21l = -xil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*xjl

          dery21i = pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*yij + t4*t5*sij*zij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat)
          dery21j = -pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*yij + t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*yjl
          dery21l = -yil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*yjl

          derz21i = pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*zij + t4*t5*sij*yij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat)
          derz21j = -pijampl*(ovrlp(iorb + 2, jorb + 1)*t5*zij + t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*zjl
          derz21l = -zil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 1)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 1)*deramplitude(jat)*zjl

          derx02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*xij + t4*t5*sij*zij) + &
                    xil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          derx02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*xij + t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*xjl
          derx02l = -xil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*xjl

          dery02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          dery02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*yjl
          dery02l = -yil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*yjl

          derz02i = pijampl*(ovrlp(iorb, jorb + 2)*t5*zij + t4*t5*sij*xij) + &
                    zil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat)
          derz02j = -pijampl*(ovrlp(iorb, jorb + 2)*t5*zij + t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*zjl
          derz02l = -zil*deramplitude(iat)*ovrlp(iorb, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb, jorb + 2)*deramplitude(jat)*zjl

          derx12i = pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat)
          derx12j = -pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*xjl
          derx12l = -xil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*xjl

          dery12i = pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*yij + t4*t5*sij*zij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat)
          dery12j = -pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*yij + t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*yjl
          dery12l = -yil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*yjl

          derz12i = pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*zij + t4*t5*sij*yij) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat)
          derz12j = -pijampl*(ovrlp(iorb + 1, jorb + 2)*t5*zij + t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*zjl
          derz12l = -zil*deramplitude(iat)*ovrlp(iorb + 1, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 1, jorb + 2)*deramplitude(jat)*zjl

          derx22i = pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*xij) + &
                    xil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat)
          derx22j = -pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*xjl
          derx22l = -xil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*xjl

          dery22i = pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*yij) + &
                    yil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat)
          dery22j = -pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*yjl
          dery22l = -yil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*yjl

          derz22i = pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*zij + t4*t5*sij*zij*2.d0) + &
                    zil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat)
          derz22j = -pijampl*(ovrlp(iorb + 2, jorb + 2)*t5*zij + t4*t5*sij*zij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*zjl
          derz22l = -zil*deramplitude(iat)*ovrlp(iorb + 2, jorb + 2)*amplitude(jat) - &
                    amplitude(iat)*ovrlp(iorb + 2, jorb + 2)*deramplitude(jat)*zjl

          do korb = 1, norb
            kkorb = norb - korb + 1  ! Also derivative in decreasing order of eigenvalues
            pevec00 = evec(korb, iorb + 0)*evec(korb, jorb + 0)
            pevec10 = evec(korb, iorb + 1)*evec(korb, jorb + 0)
            pevec20 = evec(korb, iorb + 2)*evec(korb, jorb + 0)

            pevec01 = evec(korb, iorb + 0)*evec(korb, jorb + 1)
            pevec11 = evec(korb, iorb + 1)*evec(korb, jorb + 1)
            pevec21 = evec(korb, iorb + 2)*evec(korb, jorb + 1)

            pevec02 = evec(korb, iorb + 0)*evec(korb, jorb + 2)
            pevec12 = evec(korb, iorb + 1)*evec(korb, jorb + 2)
            pevec22 = evec(korb, iorb + 2)*evec(korb, jorb + 2)

            devaldr(kkorb, 1, iat) = devaldr(kkorb, 1, iat) + pevec00*derx00i + pevec10*derx10i + pevec20*derx20i &
                                     + pevec01*derx01i + pevec11*derx11i + pevec21*derx21i &
                                     + pevec02*derx02i + pevec12*derx12i + pevec22*derx22i
            devaldr(kkorb, 1, jat) = devaldr(kkorb, 1, jat) + pevec00*derx00j + pevec10*derx10j + pevec20*derx20j &
                                     + pevec01*derx01j + pevec11*derx11j + pevec21*derx21j &
                                     + pevec02*derx02j + pevec12*derx12j + pevec22*derx22j
            devaldr(kkorb, 1, lat) = devaldr(kkorb, 1, lat) + pevec00*derx00l + pevec10*derx10l + pevec20*derx20l &
                                     + pevec01*derx01l + pevec11*derx11l + pevec21*derx21l &
                                     + pevec02*derx02l + pevec12*derx12l + pevec22*derx22l

            devaldr(kkorb, 2, iat) = devaldr(kkorb, 2, iat) + pevec00*dery00i + pevec10*dery10i + pevec20*dery20i &
                                     + pevec01*dery01i + pevec11*dery11i + pevec21*dery21i &
                                     + pevec02*dery02i + pevec12*dery12i + pevec22*dery22i
            devaldr(kkorb, 2, jat) = devaldr(kkorb, 2, jat) + pevec00*dery00j + pevec10*dery10j + pevec20*dery20j &
                                     + pevec01*dery01j + pevec11*dery11j + pevec21*dery21j &
                                     + pevec02*dery02j + pevec12*dery12j + pevec22*dery22j
            devaldr(kkorb, 2, lat) = devaldr(kkorb, 2, lat) + pevec00*dery00l + pevec10*dery10l + pevec20*dery20l &
                                     + pevec01*dery01l + pevec11*dery11l + pevec21*dery21l &
                                     + pevec02*dery02l + pevec12*dery12l + pevec22*dery22l

            devaldr(kkorb, 3, iat) = devaldr(kkorb, 3, iat) + pevec00*derz00i + pevec10*derz10i + pevec20*derz20i &
                                     + pevec01*derz01i + pevec11*derz11i + pevec21*derz21i &
                                     + pevec02*derz02i + pevec12*derz12i + pevec22*derz22i
            devaldr(kkorb, 3, jat) = devaldr(kkorb, 3, jat) + pevec00*derz00j + pevec10*derz10j + pevec20*derz20j &
                                     + pevec01*derz01j + pevec11*derz11j + pevec21*derz21j &
                                     + pevec02*derz02j + pevec12*derz12j + pevec22*derz22j
            devaldr(kkorb, 3, lat) = devaldr(kkorb, 3, lat) + pevec00*derz00l + pevec10*derz10l + pevec20*derz20l &
                                     + pevec01*derz01l + pevec11*derz11l + pevec21*derz21l &
                                     + pevec02*derz02l + pevec12*derz12l + pevec22*derz22l
          end do
        end do
      end do
    end do
  end do

1111 continue
  call rots(norb, 3*nat, devaldr(1, 1, 1), devaldr_out(1, 1, 1))
!         call copy(norb*3*nat,devaldr(1,1,1),devaldr_out(1,1,1))
  deallocate (evecn, evec, devaldr)
  deallocate (ovrlp)
  deallocate (work)
end subroutine xyz2devaldr

subroutine copy(n, a, b)
  implicit none
  integer :: n
  real*8, dimension(n) :: a
  real*8, dimension(n) :: b
  integer :: i
  do i = 1, n
    b(i) = a(i)
  end do
  return
end subroutine copy

subroutine rots(n1, n2, a, b)
  implicit none
  integer, parameter :: lot = 32
  integer :: n1
  integer :: n2
  real*8, dimension(n1, n2) :: a
  real*8, dimension(n2, n1) :: b
  integer :: jj, i, j

  do jj = 1, n2, lot
    do i = 1, n1
      do j = jj, min(n2, jj + (lot - 1))
        b(j, i) = a(i, j)
      end do
    end do
  end do
  return
end subroutine rots

subroutine rot(n1, n2, a, b)
  implicit none
  integer :: n1
  integer :: n2
  real*8, dimension(n1, n2) :: a
  real*8, dimension(n2, n1) :: b
  integer :: i, j

  do i = 1, n1
    do j = 1, n2
      b(j, i) = a(i, j)
    end do
  end do

  return
end subroutine rot

subroutine rotb(n, a, b)
  implicit none
  integer, parameter :: lot = 16
  integer :: n
  real*8, dimension(n, n) :: a
  real*8, dimension(n, n) :: b
  integer :: ii, jj, i, j

  ! loop over blocks
  do ii = 1, n, lot
    do jj = 1, n, lot
      ! loop over elements in each block
      do i = ii, min(n, ii + (lot - 1))
        do j = jj, min(n, jj + (lot - 1))
          b(j, i) = a(i, j)
        end do
      end do
    end do
  end do

  return
end subroutine rotb

subroutine atoms_sphere(width_cutoff, nex_cutoff, lat, llat, ixyzmax, nat, natx_sphere, &
                        nat_sphere, alat, rxyz, rxyz_sphere, rcov, rcov_sphere, indat, amplitude, deramplitude)
  implicit none
  real*8 :: width_cutoff
  integer :: nex_cutoff
  integer :: lat, llat, ixyzmax, nat, natx_sphere
  real*8, dimension(3, 3) :: alat
  real*8, dimension(3, nat) :: rxyz
  real*8, dimension(nat) :: rcov
  real*8, dimension(3, natx_sphere) :: rxyz_sphere
  real*8, dimension(natx_sphere) :: rcov_sphere
  real*8, dimension(natx_sphere) :: amplitude
  real*8, dimension(natx_sphere) :: deramplitude
  integer, dimension(natx_sphere) :: indat
  integer :: nat_sphere

  real*8 :: dist2, factor_cutoff
  integer :: ix, iy, iz, jat
  real*8 :: radius_cutoff, radius_cutoff2, temp, xj, yj, zj

  radius_cutoff = sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2 = radius_cutoff**2
  factor_cutoff = 1.d0/(2.d0*nex_cutoff*width_cutoff**2)
  ! write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

  nat_sphere = 0
  do jat = 1, nat
    do ix = -ixyzmax, ixyzmax
      do iy = -ixyzmax, ixyzmax
        do iz = -ixyzmax, ixyzmax
          xj = rxyz(1, jat) + ix*alat(1, 1) + iy*alat(1, 2) + iz*alat(1, 3)
          yj = rxyz(2, jat) + ix*alat(2, 1) + iy*alat(2, 2) + iz*alat(2, 3)
          zj = rxyz(3, jat) + ix*alat(3, 1) + iy*alat(3, 2) + iz*alat(3, 3)
          dist2 = (xj - rxyz(1, lat))**2 + (yj - rxyz(2, lat))**2 + (zj - rxyz(3, lat))**2
          ! write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)

          if (dist2 .le. radius_cutoff2) then
            nat_sphere = nat_sphere + 1
            if (nat_sphere .gt. natx_sphere) stop 'enlarge natx_sphere'
            !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
            temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff - 1)
            amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
            deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

            rxyz_sphere(1, nat_sphere) = xj
            rxyz_sphere(2, nat_sphere) = yj
            rxyz_sphere(3, nat_sphere) = zj
            rcov_sphere(nat_sphere) = rcov(jat)
            indat(nat_sphere) = jat
            if (dist2 .eq. 0.d0) llat = nat_sphere
          end if
        end do
      end do
    end do
  end do
end subroutine atoms_sphere

subroutine crtovrlp(nat, rxyz, alpha, cs, cp, ns, np, ovrlp)
  implicit none
  integer :: nat
  integer :: ns
  integer :: np
  real*8 :: rxyz(3, nat)
  real*8 :: ovrlp(nat*(ns + 3*np), nat*(ns + 3*np))
  real*8 :: alpha(nat), cs(10), cp(10)

  real*8 :: ai, aj
  integer :: iat, iorb, ip, is, jat, jp, js, jorb
  real*8 :: r2, sij, t1, t2, t3, t4, t5, xi, xij, xj, yi, yij, yj, zi, zij, zj

  if (ns > 10 .or. np > 10) stop 'ns > 10   .or.  np > 10  !'

  ! 1- setup the overlap matrix

  !  <s|s>
  do jat = 1, nat
    do js = 1, ns
      jorb = (jat - 1)*(ns + 3*np) + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          !!iorb=iat+(is-1)*nat
          iorb = (iat - 1)*(ns + 3*np) + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          ovrlp(iorb, jorb) = sij

        end do
      end do
    end do
  end do

  !  <pi|sj>
  do jat = 1, nat
    do js = 1, ns

      jorb = (jat - 1)*(ns + 3*np) + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np
          !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
          iorb = (iat - 1)*(ns + 3*np) + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = -2.d0*sqrt(ai)*aj/t2
          ovrlp(iorb, jorb) = t3*xij*sij
          ovrlp(iorb + 1, jorb) = t3*yij*sij
          ovrlp(iorb + 2, jorb) = t3*zij*sij

        end do
      end do
    end do
  end do

  !  <si|pj>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*(ns + 3*np) + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          !!iorb=iat+(is-1)*nat
          iorb = (iat - 1)*(ns + 3*np) + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = +2.d0*sqrt(aj)*ai/t2
          ovrlp(iorb, jorb) = t3*xij*sij
          ovrlp(iorb, jorb + 1) = t3*yij*sij
          ovrlp(iorb, jorb + 2) = t3*zij*sij

        end do
      end do
    end do
  end do

  !  <p|p>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*(ns + 3*np) + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np
          iorb = (iat - 1)*(ns + 3*np) + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t4 = 2.d0*sqrt(t1)/t2
          t5 = -2.d0*t1/t2

          ovrlp(iorb, jorb) = t4*(1.d0 + t5*xij*xij)*sij
          ovrlp(iorb + 1, jorb) = t4*(t5*yij*xij)*sij
          ovrlp(iorb + 2, jorb) = t4*(t5*zij*xij)*sij
          ovrlp(iorb, jorb + 1) = t4*(t5*xij*yij)*sij
          ovrlp(iorb + 1, jorb + 1) = t4*(1.d0 + t5*yij*yij)*sij
          ovrlp(iorb + 2, jorb + 1) = t4*(t5*zij*yij)*sij
          ovrlp(iorb, jorb + 2) = t4*(t5*xij*zij)*sij
          ovrlp(iorb + 1, jorb + 2) = t4*(t5*yij*zij)*sij
          ovrlp(iorb + 2, jorb + 2) = t4*(1.d0 + t5*zij*zij)*sij

        end do
      end do
    end do
  end do
end subroutine crtovrlp

subroutine multampspd(nat, ovrlp, amplitude, norb, ns, np, nd, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  integer :: nd
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)

  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np + 5*nd)
      jorb = (jat - 1)*(ns + 3*np + 5*nd) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np + 5*nd)
          iorb = (iat - 1)*(ns + 3*np + 5*nd) + i
          ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multampspd

subroutine multampoff(nat, ovrlp, amplitude, norb, ns, np, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)
  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np)
      jorb = (jat - 1)*(ns + 3*np) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np)
          iorb = (iat - 1)*(ns + 3*np) + i
          if (iat .eq. jat) then
            ovrla(iorb, jorb) = ovrlp(iorb, jorb)
          else
            ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
          end if
        end do
      end do
    end do
  end do

end subroutine multampoff

subroutine multamp(nat, ovrlp, amplitude, norb, ns, np, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)
  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np)
      jorb = (jat - 1)*(ns + 3*np) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np)
          iorb = (iat - 1)*(ns + 3*np) + i
          ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multamp


subroutine diagonalizeMatrix(n, aa, eval)
  implicit none

  ! Calling arguments
  integer, intent(in) :: n
  real*8, dimension(n, n), intent(in) :: aa
  real*8, dimension(n), intent(in) :: eval
  integer :: lwork
  integer :: info

  ! Local variables
  real*8, dimension(:), allocatable:: work

  lwork = 100*n
  allocate (work(lwork))
  call dsyev('v', 'l', n, aa, n, eval, work, lwork, info)
  if (info /= 0) stop ' ERROR in dsyev'
  deallocate (work)

end subroutine diagonalizeMatrix






!
! end module symmetry_penalty
