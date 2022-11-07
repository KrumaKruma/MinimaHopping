subroutine energyandforces_bazant(nat, alat0, rxyz0, etot, fxyz, deralat, stress)
!use iso_c_binding
  !   forces-edip.f
  !   -------------

  !   Version 1.0f

  !   Force and Energy Calculation with the
  !   Environment-Dependent Interatomic Potential

  !   written by Martin Z. Bazant,
  !   Department of Physics, Harvard University
  !   April - October 1997
  !   (based on forces.c by Martin Z. Bazant, June 1994)

  !   New address (2000):
  !   Professor Martin Z. Bazant
  !   Department of Mathematics 2-363B
  !   Massachusetts Institute of Technology
  !   Cambridge, MA 02139-4307

  !   E-mail:
  !   bazant@math.mit.edu

  !   translated from c to FORTRAN
  !   by Noam Bernstein, noamb@cmt.harvard.edu, December 1997
  !   modified by Max Amsler and Stefan Goedecker

  !   COPYRIGHT NOTICE
  !   ----------------
!  Copyright (C) 2001-2002 Stefan Goedecker, CEA Grenoble
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .


  !   forces-edip, copyright 1997 by Martin Z. Bazant and Harvard University.
  !   Permission is granted to use forces-edip for academic use only, at
  !   no cost. Unauthorized sale or commerical use of this software
  !   is prohibited by United States copyright law. Any publication describing
  !   research involving this software should contain the following citations,
  !   at once and in this order, to give proper credit to the theoretical work and
  !   fitting that produced EDIP and this subroutine:

  !     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
  !     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
  !     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip,
  !            Phys. Rev. B 58, 2539 (1998).

  !   This software has been extensively tested for molecular dynamics simulations
  !   on Sun, SGI and IBM architectures, but no guarantees are made.

  !   WEBSITE
  !   -------

  !   Updated versions of this software are available at the EDIP distribution site,
  !   http://pelion.eas.harvard.edu/software/EDIP/
  !   Postscript files of related papers are available at the Kaxiras group web site
  !   in the Department of Physics at Harvard University, http://pelion.eas.harvard.edu,
  !   under 'Empirical Methods'.
  !   A description of the algorithm used in this subroutine can be found
  !   in the Ph.D. Thesis of M. Z. Bazant (1997), chapter 6, on the web at
  !   http://pelion.eas.harvard.edu/~bazant/thesis/

  !   INTERFACE
  !   ---------
  !    Internally things are computed in Angstroem and eV, The input (rxyz0,alat0) is however in bohr and
  !     the output energy in Hartree, fxyz and deralat  in Hartree/Bohr

  !     nat : number of particles
  !     rxyz0 : array (3,N) of positions
  !     alat0 : lattice vecors (each column is one vector)
  !     Etot : returned energy in eV
  !     fxyz : returned forces array (3,N) in eV/Angstroms
  !     alatder  :  derivative of energy with respsct to altice vectors

  !     neighbors(p_nbrs(i)),...,neighbors(p_nbrs(i+1)) are the
  !     atoms which are "neighbors" of atom i, using a standard Verlet
  !     neighbor list. These are a global arrays, not passed, that are declared
  !     in "edip_neighbors_include.h". This way of storing atomi! positions
  !     is not unique, and will require a customized patch to the main MD program.

  !     The parameters of the potential initialized in input_EDIP_params() are global
  !     variables declared in "edip_pot_include.h".

  !   PARAMETERS
  !   ----------

  !    par_cap_A,par_cap_B,par_rh,par_a,par_sig
  !    par_lam,par_gam,par_b,par_c,par_delta
  !    par_mu,par_Qo,par_palp,par_bet,par_alp

  !    5.6714030     2.0002804     1.2085196     3.1213820     0.5774108
  !    1.4533108     1.1247945     3.1213820     2.5609104    78.7590539
  !    0.6966326   312.1341346     1.4074424     0.0070975     3.1083847

  !    Connection between these parameters and those given in the paper,
  !    Justo et al., Phys. Rev. B 58, 2539 (1998):

  !    A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))

  !    so in the paper (')
  !    A' = A*palp
  !    B' = B * palp**(-1/rh)
  !    eta = detla/Qo

  !    Non-adjustable parameters for tau(Z) from Ismail & Kaxiras, 1993,
  !    also in Bazant, Kaxiras, Justo, PRB (1997):

  !    u1 = -0.165799;
  !    u2 = 32.557;
  !    u3 = 0.286198;
  !    u4 = 0.66;
  implicit none

  !  ------------------------- VARIABLE DECLARATIONS -------------------------

  integer, intent(in) :: nat
  real*8, intent(in) :: alat0(3,3)
  real*8, intent(in) :: rxyz0(3,nat)

  real*8, intent(out) :: etot
  real*8, intent(out) :: fxyz(3,nat)
  real*8, intent(out) :: deralat(3,3)
  real*8, intent(out) :: stress(3, 3)



  integer, parameter :: nnbrx = 50 !number of geighbors
  real*8 :: Ha_eV, Bohr_Ang

  integer i, j, k, l, n, n2, n3, nz, iat
  real*8 :: dx, dy, dz, r, asqr
  real*8 :: rinv, rmainv, xinv, xinv3, den, Z, fZ
  real*8 :: dV2j, dV2ijx, dV2ijy, dV2ijz, pZ, dp
  real*8 :: temp0, temp1
  real*8 :: Qort, muhalf, u5
  real*8 :: rmbinv, winv, dwinv, tau, dtau, lcos, x, H, dHdx, dhdl
  real*8 :: dV3rij, dV3rijx, dV3rijy, dV3rijz
  real*8 :: dV3rik, dV3rikx, dV3riky, dV3rikz
  real*8 :: dV3l, dV3ljx, dV3ljy, dV3ljz, dV3lkx, dV3lky, dV3lkz
  real*8 :: dV2dZ, dxdZ, dV3dZ
  real*8 :: dEdrl, dEdrlx, dEdrly, dEdrlz
  real*8 :: bmc, cmbinv
  real*8 :: fjx, fjy, fjz, fkx, fky, fkz

  real*8 s2_t0(nnbrx), s2_t1(nnbrx), s2_t2(nnbrx), s2_t3(nnbrx), s2_dx(nnbrx), s2_dy(nnbrx), s2_dz(nnbrx), s2_r(nnbrx), &
    s3_g(nnbrx), s3_dg(nnbrx), s3_rinv(nnbrx), s3_dx(nnbrx), s3_dy(nnbrx), s3_dz(nnbrx), s3_r(nnbrx), &
    sz_df(nnbrx), sz_sum(nnbrx), sz_dx(nnbrx), sz_dy(nnbrx), sz_dz(nnbrx), sz_r(nnbrx)
  integer num2(nnbrx), num3(nnbrx), numz(nnbrx)

  integer nj, nk, nl
  !   indices for the store arrays
  real*8 ::  virial, virial_xyz(3)

  real*8 ::    par_cap_A, par_cap_B, par_rh, par_a, par_sig
  real*8 ::    par_lam, par_gam, par_b, par_c, par_delta
  real*8 ::    par_mu, par_Qo, par_palp, par_bet, par_alp
  real*8 ::    u1, u2, u3, u4
  real*8 :: par_bg
  real*8 :: par_eta
  real*8 :: cutoff
  real*8 :: delta_safe

  ! My variables
  integer :: lsta(2, nat), lstb(nnbrx*nat)
  real*8  :: rel(5, nat*nnbrx)
  real*8 :: ener_iat, ener, ener2
  real*8 :: alat(3, 3)
  real*8 :: rxyz(3, nat), alatinv(3, 3)
  real*8 :: si1_sj1, si2_sj2, si3_sj3
  real(8) :: vol

  !         End my variables
  ! EDIP parameters
  !         taken from Justo et al., Phys. Rev. B 58, 2539 (1998).
  par_cap_A = 5.6714030d0
  par_cap_B = 2.0002804d0
  par_rh = 1.2085196d0
  par_a = 3.1213820d0
  par_sig = 0.5774108d0
  par_lam = 1.4533108d0
  par_gam = 1.1247945d0
  par_b = 3.1213820d0
  par_c = 2.5609104d0
  par_delta = 78.7590539d0
  par_mu = 0.6966326d0
  par_Qo = 312.1341346d0
  par_palp = 1.4074424d0
  par_bet = 0.0070975d0
  par_alp = 3.1083847d0

  par_bg = par_a
  par_eta = par_delta/par_Qo
  cutoff = par_a  ! -1.d10
  delta_safe = 0.2d0

  u1 = -0.165799d0
  u2 = 32.557d0
  u3 = 0.286198d0
  u4 = 0.66d0
  !end parameters



  !Ha_eV = 27.211399d0
  !Bohr_Ang = 0.529177d0

  Ha_eV = 1.d0
  Bohr_Ang = 1.d0

  !Do some preparation including the construction of the pair list
  do j = 1, 3
  do i = 1, 3
    alat(i, j) = alat0(i, j)*Bohr_Ang
  end do
  end do

  do iat = 1, nat
    rxyz(1, iat) = rxyz0(1, iat)*Bohr_Ang
    rxyz(2, iat) = rxyz0(2, iat)*Bohr_Ang
    rxyz(3, iat) = rxyz0(3, iat)*Bohr_Ang
  end do
  call back2cell(nat, rxyz, alat)
  call nnlist(nat, nnbrx, alat, cutoff, rxyz, lsta, lstb, rel)
  call invertalat(alat, alatinv)

  !Allocation of temporary arrays

  !          L_x_div_2 = L_x/2.0D0
  !          L_y_div_2 = L_y/2.0D0
  !          L_z_div_2 = L_z/2.0D0

  !          do i=1, N_own
  fxyz = 0.d0
  ener = 0.d0
  ener2 = 0.d0

  virial = 0.0d0
  virial_xyz(:) = 0.0d0
  deralat = 0.d0

  !   COMBINE COEFFICIENTS

  asqr = par_a*par_a
  Qort = sqrt(par_Qo)
  muhalf = par_mu*0.5D0
  u5 = u2*u4
  bmc = par_b - par_c
  cmbinv = 1.0D0/(par_c - par_b)

  !  --- LEVEL 1: OUTER LOOP OVER ATOMS ---

  do i = 1, nat

    !   RESET COORDINATION AND NEIGHBOR NUMBERS
    ener_iat = 0.d0
    Z = 0.0d0
    n2 = 1
    n3 = 1
    nz = 1

    !  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---

    !            do n=p_nbrs(i), p_nbrs(i+1)-1
    !              j = neighbors(n)
    do n = lsta(1, i), lsta(2, i)
      j = lstb(n)

      !   PARTS OF TWO-BODY INTERACTION r<par_a

      num2(n2) = j
      !                rinv = 1.0/r
      !                dx = dx * rinv
      !                dy = dy * rinv
      !                dz = dz * rinv
      dx = -rel(1, n)
      dy = -rel(2, n)
      dz = -rel(3, n)
      r = rel(4, n)
      rinv = rel(5, n)

      rmainv = 1.d0/(r - par_a)
      s2_t0(n2) = par_cap_A*dexp(par_sig*rmainv)
      s2_t1(n2) = (par_cap_B*rinv)**par_rh
      s2_t2(n2) = par_rh*rinv
      s2_t3(n2) = par_sig*rmainv*rmainv
      s2_dx(n2) = dx
      s2_dy(n2) = dy
      s2_dz(n2) = dz
      s2_r(n2) = r
      n2 = n2 + 1

      if (n2 .gt. nnbrx) then
        stop 'WARNING1 enlarge nnbrx'
      end if

      !!Additional part from stefan
      !! coordination number calculated with soft cutoff between first and
      !! second nearest neighbor
      !        if (r.le.2.36d0) then
      !        coord_iat=coord_iat+1.d0
      !        else if (r.ge.3.83d0) then
      !        else
      !        xarg=(r-2.36d0)*(1.d0/(3.83d0-2.36d0))
      !        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
      !        endif
      !-----------------------------

      !   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b

      if (r < par_bg) then

        num3(n3) = j
        rmbinv = 1.0d0/(r - par_bg)
        temp1 = par_gam*rmbinv
        temp0 = dexp(temp1)
        s3_g(n3) = temp0
        s3_dg(n3) = -rmbinv*temp1*temp0
        s3_dx(n3) = dx
        s3_dy(n3) = dy
        s3_dz(n3) = dz
        s3_rinv(n3) = rinv
        s3_r(n3) = r
        n3 = n3 + 1

        !                  if(fixZ .eq. 0) then
        !Additional part from Stefan
        if (n3 .gt. nnbrx) then
          stop 'WARNING2 enlarge nnbrx'
        end if
        !--------------------------

        !   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b

        if (r .lt. par_b) then
          if (r .lt. par_c) then
            Z = Z + 1.0d0
          else
            xinv = bmc/(r - par_c)
            xinv3 = xinv*xinv*xinv
            den = 1.0d0/(1.d0 - xinv3)
            temp1 = par_alp*den
            fZ = dexp(temp1)
            Z = Z + fZ
            numz(nz) = j
            sz_df(nz) = fZ*temp1*den*3.0d0*xinv3*xinv*cmbinv
            !   df/dr
            sz_dx(nz) = dx
            sz_dy(nz) = dy
            sz_dz(nz) = dz
            sz_r(nz) = r
            nz = nz + 1
            !Additional part from Stefan
            if (nz .gt. nnbrx) then
              stop 'WARNING3 enlarge nnbrx'
            end if
            !---------------------
          end if
          !  r < par_C
        end if
        !  r < par_b
      end if
      !  fixZ .eq. 0
      !                end if
      !  r < par_bg
      !              end if
      !  rsqr < asqr
      !              end if
      !  dz < par_a
      !              end if
      !  dy < par_a
      !              end if
      !  dz < par_a
    end do

    !            if(fixZ .ne. 0) then
    !
    !              Z = tricks_Zfix
    !              pZ = par_palp*dexp(-par_bet*Z*Z)
    !              dp = 0.0

    !            else

    !   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES

    do nl = 1, nz - 1
      sz_sum(nl) = 0.0d0
    end do

    !   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION

    temp0 = par_bet*Z
    pZ = par_palp*dexp(-temp0*Z)
    !   bond order
    dp = -2.0d0*temp0*pZ
    !   derivative of bond order

    !            end if

    !  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---

    do nj = 1, n2 - 1

      temp0 = s2_t1(nj) - pZ

      !   two-body energy V2(rij,Z)

      ener_iat = ener_iat + temp0*s2_t0(nj)

      !   two-body forces

      dV2j = -s2_t0(nj)*(s2_t1(nj)*s2_t2(nj) + temp0*s2_t3(nj))
      !   dV2/dr
      dV2ijx = dV2j*s2_dx(nj)
      dV2ijy = dV2j*s2_dy(nj)
      dV2ijz = dV2j*s2_dz(nj)
      fxyz(1, i) = fxyz(1, i) + dV2ijx
      fxyz(2, i) = fxyz(2, i) + dV2ijy
      fxyz(3, i) = fxyz(3, i) + dV2ijz
      j = num2(nj)
      fxyz(1, j) = fxyz(1, j) - dV2ijx
      fxyz(2, j) = fxyz(2, j) - dV2ijy
      fxyz(3, j) = fxyz(3, j) - dV2ijz

      !   dV2/dr contribution to virial

      virial_xyz(1) = virial_xyz(1) - s2_r(nj)*(dV2ijx*s2_dx(nj))
      virial_xyz(2) = virial_xyz(2) - s2_r(nj)*(dV2ijy*s2_dy(nj))
      virial_xyz(3) = virial_xyz(3) - s2_r(nj)*(dV2ijz*s2_dz(nj))
      virial = virial - s2_r(nj)*(dV2ijx*s2_dx(nj) + dV2ijy*s2_dy(nj) + dV2ijz*s2_dz(nj))

      !Cell gradient part
      !My own implementation
      si1_sj1 = alatinv(1, 1)*s2_dx(nj) + alatinv(1, 2)*s2_dy(nj) + alatinv(1, 3)*s2_dz(nj)
      si2_sj2 = alatinv(2, 1)*s2_dx(nj) + alatinv(2, 2)*s2_dy(nj) + alatinv(2, 3)*s2_dz(nj)
      si3_sj3 = alatinv(3, 1)*s2_dx(nj) + alatinv(3, 2)*s2_dy(nj) + alatinv(3, 3)*s2_dz(nj)
      deralat(1, 1) = deralat(1, 1) - s2_r(nj)*dV2ijx*si1_sj1
      deralat(1, 2) = deralat(1, 2) - s2_r(nj)*dV2ijx*si2_sj2
      deralat(1, 3) = deralat(1, 3) - s2_r(nj)*dV2ijx*si3_sj3
      deralat(2, 1) = deralat(2, 1) - s2_r(nj)*dV2ijy*si1_sj1
      deralat(2, 2) = deralat(2, 2) - s2_r(nj)*dV2ijy*si2_sj2
      deralat(2, 3) = deralat(2, 3) - s2_r(nj)*dV2ijy*si3_sj3
      deralat(3, 1) = deralat(3, 1) - s2_r(nj)*dV2ijz*si1_sj1
      deralat(3, 2) = deralat(3, 2) - s2_r(nj)*dV2ijz*si2_sj2
      deralat(3, 3) = deralat(3, 3) - s2_r(nj)*dV2ijz*si3_sj3

      !              if(fixZ .eq. 0) then

      !  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---

      dV2dZ = -dp*s2_t0(nj)
      do nl = 1, nz - 1
        sz_sum(nl) = sz_sum(nl) + dV2dZ
      end do

      !              end if
      !  fixZ
    end do
    !Commented out by Stefan
    !            if(fixZ .ne. 0) then
    !              winv = Qort*dexp(-muhalf*Z)
    !              dwinv = 0.0
    !              temp0 = dexp(-u4*Z)
    !              tau = u1+u2*temp0*(u3-temp0)
    !              dtau = 0.0
    !            else
    !----------------------------------

    !   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION

    winv = Qort*exp(-muhalf*Z)
    !   inverse width of angular function
    dwinv = -muhalf*winv
    !   its derivative
    temp0 = exp(-u4*Z)
    tau = u1 + u2*temp0*(u3 - temp0)
    !   -cosine of angular minimum
    dtau = u5*temp0*(2.d0*temp0 - u3)
    !   its derivative
    !            end if

    !  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---

    do nj = 1, n3 - 2

      j = num3(nj)

      !  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---

      do nk = nj + 1, n3 - 1

        k = num3(nk)

        !   angular function h(l,Z)

        lcos = s3_dx(nj)*s3_dx(nk) + s3_dy(nj)*s3_dy(nk) + s3_dz(nj)*s3_dz(nk)
        x = (lcos + tau)*winv
        temp0 = exp(-x*x)

        H = par_lam*(1.d0 - temp0 + par_eta*x*x)
        dHdx = 2.d0*par_lam*x*(temp0 + par_eta)

        dhdl = dHdx*winv

        !   three-body energy

        temp1 = s3_g(nj)*s3_g(nk)
        ener_iat = ener_iat + temp1*H

        !   (-) radial force on atom j

        dV3rij = s3_dg(nj)*s3_g(nk)*H
        dV3rijx = dV3rij*s3_dx(nj)
        dV3rijy = dV3rij*s3_dy(nj)
        dV3rijz = dV3rij*s3_dz(nj)
        fjx = dV3rijx
        fjy = dV3rijy
        fjz = dV3rijz

        !   (-) radial force on atom k

        dV3rik = s3_g(nj)*s3_dg(nk)*H
        dV3rikx = dV3rik*s3_dx(nk)
        dV3riky = dV3rik*s3_dy(nk)
        dV3rikz = dV3rik*s3_dz(nk)
        fkx = dV3rikx
        fky = dV3riky
        fkz = dV3rikz

        !   (-) angular force on j

        dV3l = temp1*dhdl
        dV3ljx = dV3l*(s3_dx(nk) - lcos*s3_dx(nj))*s3_rinv(nj)
        dV3ljy = dV3l*(s3_dy(nk) - lcos*s3_dy(nj))*s3_rinv(nj)
        dV3ljz = dV3l*(s3_dz(nk) - lcos*s3_dz(nj))*s3_rinv(nj)
        fjx = fjx + dV3ljx
        fjy = fjy + dV3ljy
        fjz = fjz + dV3ljz

        !   (-) angular force on k

        dV3lkx = dV3l*(s3_dx(nj) - lcos*s3_dx(nk))*s3_rinv(nk)
        dV3lky = dV3l*(s3_dy(nj) - lcos*s3_dy(nk))*s3_rinv(nk)
        dV3lkz = dV3l*(s3_dz(nj) - lcos*s3_dz(nk))*s3_rinv(nk)
        fkx = fkx + dV3lkx
        fky = fky + dV3lky
        fkz = fkz + dV3lkz

        !   apply radial + angular forces to i, j, k

        fxyz(1, j) = fxyz(1, j) - fjx
        fxyz(2, j) = fxyz(2, j) - fjy
        fxyz(3, j) = fxyz(3, j) - fjz
        fxyz(1, k) = fxyz(1, k) - fkx
        fxyz(2, k) = fxyz(2, k) - fky
        fxyz(3, k) = fxyz(3, k) - fkz
        fxyz(1, i) = fxyz(1, i) + fjx + fkx
        fxyz(2, i) = fxyz(2, i) + fjy + fky
        fxyz(3, i) = fxyz(3, i) + fjz + fkz

        !   dV3/dR contributions to virial

        virial = virial - s3_r(nj)*(fjx*s3_dx(nj) + fjy*s3_dy(nj) + fjz*s3_dz(nj))
        virial = virial - s3_r(nk)*(fkx*s3_dx(nk) + fky*s3_dy(nk) + fkz*s3_dz(nk))
        virial_xyz(1) = virial_xyz(1) - s3_r(nj)*(fjx*s3_dx(nj))
        virial_xyz(2) = virial_xyz(2) - s3_r(nj)*(fjy*s3_dy(nj))
        virial_xyz(3) = virial_xyz(3) - s3_r(nj)*(fjz*s3_dz(nj))
        virial_xyz(1) = virial_xyz(1) - s3_r(nk)*(fkx*s3_dx(nk))
        virial_xyz(2) = virial_xyz(2) - s3_r(nk)*(fky*s3_dy(nk))
        virial_xyz(3) = virial_xyz(3) - s3_r(nk)*(fkz*s3_dz(nk))

        !Cell gradient part
        !My own implementation
        si1_sj1 = alatinv(1, 1)*s3_dx(nj) + alatinv(1, 2)*s3_dy(nj) + alatinv(1, 3)*s3_dz(nj)
        si2_sj2 = alatinv(2, 1)*s3_dx(nj) + alatinv(2, 2)*s3_dy(nj) + alatinv(2, 3)*s3_dz(nj)
        si3_sj3 = alatinv(3, 1)*s3_dx(nj) + alatinv(3, 2)*s3_dy(nj) + alatinv(3, 3)*s3_dz(nj)
        deralat(1, 1) = deralat(1, 1) - s3_r(nj)*fjx*si1_sj1
        deralat(1, 2) = deralat(1, 2) - s3_r(nj)*fjx*si2_sj2
        deralat(1, 3) = deralat(1, 3) - s3_r(nj)*fjx*si3_sj3
        deralat(2, 1) = deralat(2, 1) - s3_r(nj)*fjy*si1_sj1
        deralat(2, 2) = deralat(2, 2) - s3_r(nj)*fjy*si2_sj2
        deralat(2, 3) = deralat(2, 3) - s3_r(nj)*fjy*si3_sj3
        deralat(3, 1) = deralat(3, 1) - s3_r(nj)*fjz*si1_sj1
        deralat(3, 2) = deralat(3, 2) - s3_r(nj)*fjz*si2_sj2
        deralat(3, 3) = deralat(3, 3) - s3_r(nj)*fjz*si3_sj3

        !Cell gradient part
        !My own implementation
        si1_sj1 = alatinv(1, 1)*s3_dx(nk) + alatinv(1, 2)*s3_dy(nk) + alatinv(1, 3)*s3_dz(nk)
        si2_sj2 = alatinv(2, 1)*s3_dx(nk) + alatinv(2, 2)*s3_dy(nk) + alatinv(2, 3)*s3_dz(nk)
        si3_sj3 = alatinv(3, 1)*s3_dx(nk) + alatinv(3, 2)*s3_dy(nk) + alatinv(3, 3)*s3_dz(nk)
        deralat(1, 1) = deralat(1, 1) - s3_r(nk)*fkx*si1_sj1
        deralat(1, 2) = deralat(1, 2) - s3_r(nk)*fkx*si2_sj2
        deralat(1, 3) = deralat(1, 3) - s3_r(nk)*fkx*si3_sj3
        deralat(2, 1) = deralat(2, 1) - s3_r(nk)*fky*si1_sj1
        deralat(2, 2) = deralat(2, 2) - s3_r(nk)*fky*si2_sj2
        deralat(2, 3) = deralat(2, 3) - s3_r(nk)*fky*si3_sj3
        deralat(3, 1) = deralat(3, 1) - s3_r(nk)*fkz*si1_sj1
        deralat(3, 2) = deralat(3, 2) - s3_r(nk)*fkz*si2_sj2
        deralat(3, 3) = deralat(3, 3) - s3_r(nk)*fkz*si3_sj3

        !                if(fixZ .eq. 0) then

        !   prefactor for 4-body forces from coordination
        dxdZ = dwinv*(lcos + tau) + winv*dtau
        dV3dZ = temp1*dHdx*dxdZ

        !  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---

        do nl = 1, nz - 1
          sz_sum(nl) = sz_sum(nl) + dV3dZ
        end do
          !                end if
      end do
    end do

    !            if(fixZ .eq. 0) then

    !  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---

    do nl = 1, nz - 1

      dEdrl = sz_sum(nl)*sz_df(nl)
      dEdrlx = dEdrl*sz_dx(nl)
      dEdrly = dEdrl*sz_dy(nl)
      dEdrlz = dEdrl*sz_dz(nl)
      fxyz(1, i) = fxyz(1, i) + dEdrlx
      fxyz(2, i) = fxyz(2, i) + dEdrly
      fxyz(3, i) = fxyz(3, i) + dEdrlz
      l = numz(nl)
      fxyz(1, l) = fxyz(1, l) - dEdrlx
      fxyz(2, l) = fxyz(2, l) - dEdrly
      fxyz(3, l) = fxyz(3, l) - dEdrlz

      !   dE/dZ*dZ/dr contribution to virial

      virial = virial - sz_r(nl)*(dEdrlx*sz_dx(nl) + dEdrly*sz_dy(nl) + dEdrlz*sz_dz(nl))
      virial_xyz(1) = virial_xyz(1) - sz_r(nl)*(dEdrlx*sz_dx(nl))
      virial_xyz(2) = virial_xyz(2) - sz_r(nl)*(dEdrly*sz_dy(nl))
      virial_xyz(3) = virial_xyz(3) - sz_r(nl)*(dEdrlz*sz_dz(nl))

      !Cell gradient part
      !My own implementation
      si1_sj1 = alatinv(1, 1)*sz_dx(nl) + alatinv(1, 2)*sz_dy(nl) + alatinv(1, 3)*sz_dz(nl)
      si2_sj2 = alatinv(2, 1)*sz_dx(nl) + alatinv(2, 2)*sz_dy(nl) + alatinv(2, 3)*sz_dz(nl)
      si3_sj3 = alatinv(3, 1)*sz_dx(nl) + alatinv(3, 2)*sz_dy(nl) + alatinv(3, 3)*sz_dz(nl)
      deralat(1, 1) = deralat(1, 1) - sz_r(nl)*dEdrlx*si1_sj1
      deralat(1, 2) = deralat(1, 2) - sz_r(nl)*dEdrlx*si2_sj2
      deralat(1, 3) = deralat(1, 3) - sz_r(nl)*dEdrlx*si3_sj3
      deralat(2, 1) = deralat(2, 1) - sz_r(nl)*dEdrly*si1_sj1
      deralat(2, 2) = deralat(2, 2) - sz_r(nl)*dEdrly*si2_sj2
      deralat(2, 3) = deralat(2, 3) - sz_r(nl)*dEdrly*si3_sj3
      deralat(3, 1) = deralat(3, 1) - sz_r(nl)*dEdrlz*si1_sj1
      deralat(3, 2) = deralat(3, 2) - sz_r(nl)*dEdrlz*si2_sj2
      deralat(3, 3) = deralat(3, 3) - sz_r(nl)*dEdrlz*si3_sj3

    end do

    !           end if
    ener = ener + ener_iat
    ener2 = ener2 + ener_iat**2
  end do

  !        call getvol(alat,vol)
  !          if(vol.lt.0.d0) then
  !            write(77,*) nat
  !            write(77,*) alat(:,1)
  !            write(77,*) alat(:,2)
  !            write(77,*) alat(:,3)
  !            do i=1,nat
  !              write(77,*)  rxyz(:,i),"Si"
  !            enddo
  !
  !          endif

  etot = ener

  etot = etot/Ha_eV
  do iat = 1, nat
    fxyz(1, iat) = -fxyz(1, iat)/Ha_eV*Bohr_Ang
    fxyz(2, iat) = -fxyz(2, iat)/Ha_eV*Bohr_Ang
    fxyz(3, iat) = -fxyz(3, iat)/Ha_eV*Bohr_Ang
  end do
  do j = 1, 3
  do i = 1, 3
    deralat(i, j) = deralat(i, j)/Ha_eV*Bohr_Ang
  end do
  end do
  ! Stress
  ! Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
  !        do i=1,3
  !           tmplat(:,i)=latvec(i,:)
  !        enddo


  !vol = (alat0(1, 1)*alat0(2, 2)*alat0(3, 3) - alat0(1, 1)*alat0(2, 3)*alat0(3, 2) - &
  !  alat0(1, 2)*alat0(2, 1)*alat0(3, 3) + alat0(1, 2)*alat0(2, 3)*alat0(3, 1) + &
  !  alat0(1, 3)*alat0(2, 1)*alat0(3, 2) - alat0(1, 3)*alat0(2, 2)*alat0(3, 1))
  !stress=-matmul(deralat,transpose(alat0))/vol

  vol = (alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
    alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
    alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))

  !print*, "F:  ", alat
  stress=-matmul(deralat,transpose(alat))/vol
  !stress = -matmul(transpose(deralat), transpose(alat))/vol
  !stress = deralat

  !stress = stress / Ha_eV *Bohr_Ang**3
  !        strten(1) = stress(1,1)
  !        strten(2) = stress(2,2)
  !        strten(3) = stress(3,3)
  !        strten(6) = stress(2,1)
  !        strten(5) = stress(3,1)
  !        strten(4) = stress(3,2)
  !!This is not very clear yet...
  !        strten=strten/Ha_eV*Bohr_Ang**3
end subroutine energyandforces_bazant

subroutine nnlist(nat, nnbrx, alat, cutoff, rxyz, lsta, lstb, rel)
  implicit none
  integer, parameter :: nwork=1000
  real(8) :: rxyz(3, nat), rel(5, nat*nnbrx)
  real(8) :: alat(3, 3)
  integer :: lsta(2, nat), lstb(nnbrx*nat)
  real(8) :: alatalat(3, 3), eigalat(3), workalat(nwork)
  integer :: ixyzmax, info, ind, iat, jat, ix, iy, iz, nat, nnbrx, i, j
  real(8) :: cutoff2, cutoff, dist2, relx, rely, relz, tt, ttinv, xj, yj, zj

  if (alat(1, 1)*alat(2, 2)*alat(3, 3) .eq. 0.d0) then ! no periodic boundary condition
    ixyzmax = 0
  else  ! periodic boundary conditions
    do i = 1, 3
    do j = 1, 3
      alatalat(i, j) = alat(1, i)*alat(1, j) + alat(2, i)*alat(2, j) + alat(3, i)*alat(3, j)
    end do
    end do
    call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
    !   write(*,*) !  'alat !  EVals',eigalat !  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
    ! ixyzmax determines over how many periodiv images one has to search to fill the sphere with atoms
    ixyzmax = int(sqrt(1.d0/eigalat(1))*cutoff) + 1
  end if
  cutoff2 = cutoff**2

  !  write(*,*) 'ixyzmax ',ixyzmax

  ind = 0
  do iat = 1, nat
    lsta(1, iat) = ind + 1

    do jat = 1, nat
      do ix = -ixyzmax, ixyzmax
        do iy = -ixyzmax, ixyzmax
          do iz = -ixyzmax, ixyzmax
            xj = rxyz(1, jat) + ix*alat(1, 1) + iy*alat(1, 2) + iz*alat(1, 3)
            yj = rxyz(2, jat) + ix*alat(2, 1) + iy*alat(2, 2) + iz*alat(2, 3)
            zj = rxyz(3, jat) + ix*alat(3, 1) + iy*alat(3, 2) + iz*alat(3, 3)
            relx = xj - rxyz(1, iat)
            rely = yj - rxyz(2, iat)
            relz = zj - rxyz(3, iat)
            dist2 = relx**2 + rely**2 + relz**2

            if (dist2 .gt. 1.d-20 .and. dist2 .le. cutoff2) then
              ind = ind + 1
              if (ind .gt. nnbrx*nat) stop 'enlarge nnbrx'
              lstb(ind) = jat
              tt = sqrt(dist2)
              ttinv = 1.d0/tt
              rel(1, ind) = relx*ttinv
              rel(2, ind) = rely*ttinv
              rel(3, ind) = relz*ttinv
              rel(4, ind) = tt
              rel(5, ind) = ttinv
            end if
          end do
        end do
      end do
    end do
    lsta(2, iat) = ind
  end do

  !  do iat=1,nat
  !  write(*,'(i3,1x,20(1x,i2))') iat,(lstb(j),j=lsta(1,iat),lsta(2,iat))
  !  write(*,'(i3,1x,20(1x,e9.2))') iat,(rel(4,j),j=lsta(1,iat),lsta(2,iat))
  !  enddo
end subroutine nnlist

subroutine invertalat(alat, alatinv)
  !Invert alat matrix
  implicit none
  real(8) :: alat(3, 3), alatinv(3, 3), div

  div = (alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
         alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
         alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))
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
end subroutine invertalat

subroutine cart2frac(nat,alat,rxyz,xyzred)
  !! converts cartesian coordinates rxyz to reduced coordinates xyzred
  implicit none
  integer, intent(in) :: nat
  !! Number of Atoms
  real*8, intent(in), dimension(3,3) :: alat
  !! Lattice Vectors.
  real*8, intent(in), dimension(3,nat) :: rxyz
  !! Position of the Atoms in cartesian coorinates.
  real*8, intent(out), dimension(3,nat) :: xyzred
  !! Position of the Atoms in reduced coordinates.
  real(8) :: alatinv(3,3), div
  integer :: iat

    div=alat(1,1)*alat(2,2)*alat(3,3)-alat(1,1)*alat(2,3)*alat(3,2)- &
        alat(1,2)*alat(2,1)*alat(3,3)+alat(1,2)*alat(2,3)*alat(3,1)+ &
        alat(1,3)*alat(2,1)*alat(3,2)-alat(1,3)*alat(2,2)*alat(3,1)
    div=1.d0/div
      alatinv(1,1) = (alat(2,2)*alat(3,3)-alat(2,3)*alat(3,2))*div
      alatinv(1,2) =-(alat(1,2)*alat(3,3)-alat(1,3)*alat(3,2))*div
      alatinv(1,3) = (alat(1,2)*alat(2,3)-alat(1,3)*alat(2,2))*div
      alatinv(2,1) =-(alat(2,1)*alat(3,3)-alat(2,3)*alat(3,1))*div
      alatinv(2,2) = (alat(1,1)*alat(3,3)-alat(1,3)*alat(3,1))*div
      alatinv(2,3) =-(alat(1,1)*alat(2,3)-alat(1,3)*alat(2,1))*div
      alatinv(3,1) = (alat(2,1)*alat(3,2)-alat(2,2)*alat(3,1))*div
      alatinv(3,2) =-(alat(1,1)*alat(3,2)-alat(1,2)*alat(3,1))*div
      alatinv(3,3) = (alat(1,1)*alat(2,2)-alat(1,2)*alat(2,1))*div

      do iat=1,nat
      xyzred(1,iat)=alatinv(1,1)*rxyz(1,iat)+alatinv(1,2)*rxyz(2,iat)+alatinv(1,3)*rxyz(3,iat)
      xyzred(2,iat)=alatinv(2,1)*rxyz(1,iat)+alatinv(2,2)*rxyz(2,iat)+alatinv(2,3)*rxyz(3,iat)
      xyzred(3,iat)=alatinv(3,1)*rxyz(1,iat)+alatinv(3,2)*rxyz(2,iat)+alatinv(3,3)*rxyz(3,iat)
      enddo
end subroutine cart2frac

subroutine frac2cart(nat, alat, xyzred, rxyz)
  !! Converts reduced coordinates xyzred to cartesian coordinates rxyz
  implicit none
  integer, intent(in) :: nat
  !! Number of atoms.
  real*8, intent(in), dimension(3,3) :: alat
  !! Lattice Vecors
  real*8, dimension(3,nat), intent(in) :: xyzred
  !! Position of the atoms in reduced coordinates.
  real*8, dimension(3,nat), intent(out) :: rxyz
  !! Position of the atoms in cartesian coordinates.
  integer :: iat, i, j
  real(8) :: t

  do iat=1,nat
    do i = 1, 3
      t = 0.d0
      do j = 1, 3
        t = t + xyzred(j,iat) * alat(i, j)
      end do
      rxyz(i,iat) = t
    end do
  enddo
end subroutine frac2cart

subroutine back2cell(nat,rxyz,alat)
  !! Translates atoms outside the cell back into the cell.
  implicit none
  integer, intent(in) :: nat
  !! Number of atoms
  real*8, dimension(3,nat), intent(inout) :: rxyz
  !! Positions of the atoms.
  real*8, dimension(3,3), intent(in) :: alat
  !! Lattice vectors of the atoms.
  real(8) :: xyzred(3,nat)
  integer :: iat, l

  call cart2frac(nat,alat,rxyz,xyzred)
  do iat=1,nat
    do l=1,3
      xyzred(l,iat)=modulo(xyzred(l,iat),1.d0)
    enddo
  enddo
  call frac2cart(nat,alat,xyzred,rxyz)
end subroutine back2cell
