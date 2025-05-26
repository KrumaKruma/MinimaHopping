    implicit real*8 (a-h,o-z)
    parameter (nat=12)
    dimension rxyz(3,nat), rcov(nat)
    Character*2 atomname(nat)

    open(unit=11,file="posinp.xyz")
    read(11,*) natp
    read(11, *)
    if (nat.ne.natp) stop 'nat'
    do iat=1,nat
    read(11,*) atomname(iat),(rxyz(j,iat),j=1,3)
    if (atomname(iat).eq.'Si') rcov(iat)=1.11d0
    if (atomname(iat).eq.'O') rcov(iat)=0.66d0
!    write(*,*)  atomname(iat),(rxyz(j,iat),j=1,3), rcov(iat)
    !  rxyz(:,iat)= 2.*rxyz(:,iat)
    enddo

    open(unit=22,file="frags.xyz")

    call fixfrag(nat,rxyz,rcov,atomname)


    end


    subroutine fixfrag(nat,rxyz,rcov,atomname)
! nat: number of atoms
! checks whether a cluster is fragmented (distance of the closest atoms more than 1.3 times the sum of the covalent radii)
! rxyz: possibly fragmented position  on input, fixed position on output
! rcov: covalent radius of the atoms (in the same units as rxyz)
! atomname: array with the chemical symbols (needed only for visualization of the steps done by the subroutine)
    implicit real*8 (a-h,o-z)
    dimension rxyz(3,nat), rcov(nat), vec(3)
    Character*2 atomname(nat)
    ! belong is 'true' if the atom belongs to fragment containing the first atom and 'false' if it belongs to other fragments
    logical belong(nat),logexit,debug ,loggrow

    debug=.true.
    icount=0

 do ! iterate in case there are several fragments; bring back one fragment in each iteration
     if (debug) then 
     write(*,*) 'next iteration'
     write(22,*) nat
     write(22,*)  icount
     do iat=1,nat
     write(22,*) atomname(iat),(rxyz(j,iat),j=1,3)
     enddo
     endif
    icount=icount+1

    belong(1)=.true.   ! to start with, the first fragment just contains atom 1
    do iat=2,nat   
    belong(iat)=.false.
    enddo
100 continue   

! find atoms belonging to the first fragment
    loggrow=.false.
    do  iiat=1,nat
        if (belong(iiat) .eqv. .true.) then
          do iat=1,iiat-1
          d2=(rxyz(1,iat)-rxyz(1,iiat))**2 + (rxyz(2,iat)-rxyz(2,iiat))**2 + (rxyz(3,iat)-rxyz(3,iiat))**2  
          if (d2.le.((rcov(iiat)+rcov(iat))*1.3d0)**2)  then 
              if (belong(iat) .eqv. .false.) loggrow=.true.   ! new atom has been assigned to fragment
              belong(iat)=.true.
          endif
          enddo
          do iat=iiat+1,nat
          d2=(rxyz(1,iat)-rxyz(1,iiat))**2 + (rxyz(2,iat)-rxyz(2,iiat))**2 + (rxyz(3,iat)-rxyz(3,iiat))**2  
          if (d2.le.((rcov(iiat)+rcov(iat))*1.3d0)**2)  then 
              if (belong(iat) .eqv. .false.) loggrow=.true.   ! new atom has been assigned to fragment
              belong(iat)=.true.
          endif
          enddo
        endif
    enddo

       if (debug) then
       write(*,*) 'updatbel',belong
       write(*,*) ' '
       endif
          if (loggrow) goto 100  ! keep growing the fragment because last growth step still gave a larger fragment

       if (debug) then
       write(*,*) ' belong ',belong
       endif


     logexit = .true.
     do iat=2,nat
     if (belong(iat) .eqv. .false.) logexit = .false. ! some more fragments have to be fixed; continue iterations
     enddo
     if (debug) write(*,*) 'logexit ',logexit
     if (logexit) goto 1000

! All the atoms in the first fragment have been found
! Now find shortest distance between atoms in the two fragments
     dmin=1.d100
     do iat=1,nat
     if (belong(iat)) then
     do jat=1,nat
     if (belong(jat) .eqv. .false.) then
      d2=(rxyz(1,iat)-rxyz(1,jat))**2 + (rxyz(2,iat)-rxyz(2,jat))**2 + (rxyz(3,iat)-rxyz(3,jat))**2  
      if (d2.lt.dmin) then
        dmin=d2
        iiat=iat
        jjat=jat
      endif
     endif
     enddo
     endif
     enddo

     if (debug) Write(*,*) 'shortest distance ',iiat,jjat,sqrt(dmin)
    vec(1)=rxyz(1,iiat)-rxyz(1,jjat)
    vec(2)=rxyz(2,iiat)-rxyz(2,jjat)
    vec(3)=rxyz(3,iiat)-rxyz(3,jjat)

    ! find a translation vec such that the covalent radii of the closest atoms touch
    veclength=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
    scale=(veclength-rcov(iiat)-rcov(jjat))/veclength
    vec(1)=vec(1)*scale
    vec(2)=vec(2)*scale
    vec(3)=vec(3)*scale

    do iat=1,nat
    if (belong(iat).eqv. .false.) then
        rxyz(1,iat)=rxyz(1,iat)+vec(1)
        rxyz(2,iat)=rxyz(2,iat)+vec(2)
        rxyz(3,iat)=rxyz(3,iat)+vec(3)
    endif
    enddo

 enddo
1000 continue

! write fixed positions
    if (debug) then
    write(22,*) nat
    write(22,*)  icount
    do iat=1,nat
    write(22,*) atomname(iat),(rxyz(j,iat),j=1,3)
    enddo
    endif


    end subroutine
