program mixdyn
! **********************************************************************
!
! TIME INTEGRATION IMPLICIT-EXPLICITY ALGORITHM
!
! Variables
!
!   NPOIN   Number of nodes
!   NELEM   Number of elements
!   NDOFN   Number of degrees of freedom per node
!   NMATS   Number of different material sets
!   NVFIX   Number of bound nodes
!   NTYPE   Type of problem  1: plane stress, 2: plane strain, 3: axisymmetric problem
!   NNODE   Number of nodes for element (4, 8, 9)    
!   NPROP   Number of material properties (11)
!   NGAUS   Integration rule for stiffnes matrix
!   NDIME   Number of coordinates (2)    
!   NSTRE   Number of stress component (3 plain strees strain, 4 axisymmetric))   
!   NCRIT   Yield criterion (1=Tresca, 2=VonMises, 3=Mohr-Coulomb, 4=Drucker-Prager)
!   NPREV           
!   NCONM   Number of concentrated masses
!   NLAPS   Analysis type indicator (0=elastic analysis, 
!                                    1=elasto-plastic small displacement analysis, 
!                                    2=elastic large displacement analysis )
!   NGAUM   Integration rule for mass matrix 
!   NRADS   (r, z) coordinates, (R,Theta) coordinates
!   LNODS   Nodal connection numbers    
!   COORD   Nodal coordinate
!   IFPRE   Nodal constraints
!   PROPS   Material properties     PROPS(NUMAT, 1)= Young's modulus
!                                   PROPS(NUMAT, 2)= Poisson's ratio
!                                   PROPS(NUMAT, 3)= thickness for plane stress
!                                   PROPS(NUMAT, 4)= mass density
!                                   PROPS(NUMAT, 5)= temperature coefficient
!                                   PROPS(NUMAT, 6)= reference yield value
!                                   PROPS(NUMAT, 7)= hardening parameter
!                                   PROPS(NUMAT, 8)= friction angle
!                                   PROPS(NUMAT, 9)= fluidity parameter
!                                   PROPS(NUMAT,10)= exponent
!                                   PROPS(NUMAT,11)= NFLOW code (1=power law, <> 1 exponential law)     
!   NSTEP   Total number of time steps
!   NOUTD   Number of steps for printing on tape 10 and 11 of the selected displacement
!   NOUTP 
!   NREQD   Number of nodes for selective output of displacements at NOUTD steps.
!   NREQS   Number of integration points for selective output of stresses at every NOUTP step.
!   NACCE   Number of acceleration ordinates(if IFUNC # 0, NACCE is not used)
!   IFUNC   Time function code. IFUNC=0 Acceleration time history
!                               IFUNC=1 Heavside function
!                               IFUNC=2 Harmonic exicitation
!   IFIXD   Indicator for exicitation.  IFIXD=0 : horizontal acceleration read fron tape 7,
!                                                 vertical acceleration read fron tape 12
!                                       IFIXD=1 : vertical acceleration read fron tape 12
!                                       IFIXD=2 : horizontal acceleration read fron tape 7
!   MITER   Maximum number of iterations
!   KSTEP   Number of steps after which the stiffness matrix is reformed
!   IPRED   = 1 standard algorithm, = 2 modified algorithm
!
!   DTIME   Time step length
!   DTEND   Time at the end of the excitation force.
!   DTREC   Time step ofacceleration records
!   AALFA   Damping parameter (C=aalfa*M, aalfa=2*csi*omega)
!   BEETA   Damping parameter (C=beeta*M)
!   DELTA   Newmark's integration parameter (delta=0.25*(gamma+ 0.5)^2)
!   GAAMA   Newmark's integration parameter (gamma > 0.5 for stable solution)    
!   AZERO   Constant for armonic exicitatoin    
!   BZERO   Constant for arminic exicitation
!   OMEGA   Constant for armonic exicitation
!   TOLER   Tollerance
!   NPRQD   Nodal points at which displacement history is required
!   NGRQS   Integration points at which stress history is required
!   INTGR   Implicit-explicit element indicator
!   DISPI   Initial nodal displacements    
!   VELOI   initial nodal velocity
!   FORCE   Previous load state         TODO to check        
!   STRIN   Previous stress state       TODO to check    
!   IPLOD   Point load indicator
!   IGRAV   Gravity load indicator
!   IEDGE   Edge load indicator
!   ITEMP   Temperature load indicator
!   THETA   Angle of gravity axis to the positive y axis
!   GRAVY   Gravity constant
!   NEDGE   Number of loaded edge
!   NEASS   Element number with edge load
!   NOPRS   Edge nodes read in anticlockwise sequence
!   PRESS   PRESS(*, 1) Normal component of edge load, PRESS(*, 2) Tangential component of edge load for each node.
!   NODPT   Node number
!   TEMPE   Nodal temperature
!   IPOIN   Current nodal point with concentrated mass.            
!   XCMAS   Concentrated mass associated with the x-direction
!   YCMAS   Concentrated mass associated with the y-direction 
!       
! **********************************************************************
    use model

    implicit none

    integer :: istep, iiter, nchek

    real :: consd, consf

    integer :: dt(8)

    ! Start
    call date_and_time(values=dt)
    write(*, '(a7, i4, 5(a, i2.2), a, i3/)') 'Start ', dt(1), '-', dt(2), '-', dt(3), ' ', dt(5), ':', dt(6), ':', dt(7), '.', dt(8)

    ! Open the file I/O
    open(5, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/input.dat',  status='OLD')
    open(6, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/result.dat', status='OLD')

    call date_and_time(values=dt)
    write(6, '(a7, i4, 5(a, i2.2), a, i3/)') 'Start ', dt(1), '-', dt(2), '-', dt(3), ' ', dt(5), ':', dt(6), ':', dt(7), '.', dt(8)

    call contol()
    call inputd()
    call intime()
    goto 100
    call prevos()
    call loadpl()
    call lumass()
    call linkin()
    do istep = 1, nstep
        do iiter=1, miter
            call gstiff(istep)
            call impexp(consd, consf, iiter, istep)
            call resepl(rload, iiter, resid, istep)
            call itrate(consd, consf, iiter)
            if(nchek .eq. 1) exit
        end do
        call outdyn(iiter, istep)
    end do

!   Close file I/O
    close(5)
    close(6)

100 continue
    stop

end program mixdyn
