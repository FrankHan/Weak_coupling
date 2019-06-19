!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Chen Zhang, Michigan State University
!> @brief crystallite state integration functions and reporting of results
!--------------------------------------------------------------------------------------------------

module crystallite
  use prec, only: &
    pReal, &
    pStringLen
  use rotations, only: &
    rotation
  use FEsolving, only:  &
    FEsolving_execElem, &
    FEsolving_execIP
  use material, only: &
    homogenization_Ngrains
  use future
 
  implicit none
 
  private
  character(len=64),         dimension(:,:),          allocatable, private :: &
    crystallite_output                                                                              !< name of each post result output
  integer,                                                         public, protected :: &
    crystallite_maxSizePostResults                                                                  !< description not available
  integer,                   dimension(:),            allocatable, public, protected :: &
    crystallite_sizePostResults                                                                     !< description not available
  integer,                   dimension(:,:),          allocatable, private :: &
    crystallite_sizePostResult                                                                      !< description not available
 
  real(pReal),               dimension(:,:,:),        allocatable, public :: &
    crystallite_dt                                                                                  !< requested time increment of each grain
  real(pReal),               dimension(:,:,:),        allocatable, private :: &
    crystallite_subdt, &                                                                            !< substepped time increment of each grain
    crystallite_subFrac, &                                                                          !< already calculated fraction of increment
    crystallite_subStep                                                                             !< size of next integration step
  type(rotation),            dimension(:,:,:),        allocatable, private :: &
    crystallite_orientation, &                                                                      !< orientation 
    crystallite_orientation0                                                                        !< initial orientation
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public, protected :: &
    crystallite_Fe, &                                                                               !< current "elastic" def grad (end of converged time step)
    crystallite_P                                                                                   !< 1st Piola-Kirchhoff stress per grain
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_S, &                                                                                !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_S0, &                                                                               !< 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partionedS0, &                                                                      !< 2nd Piola-Kirchhoff stress vector at start of homog inc
    crystallite_Fp, &                                                                               !< current plastic def grad (end of converged time step)
    crystallite_Fp0, &                                                                              !< plastic def grad at start of FE inc
    crystallite_partionedFp0,&                                                                      !< plastic def grad at start of homog inc
    crystallite_Fi, &                                                                               !< current intermediate def grad (end of converged time step)
    crystallite_Fi0, &                                                                              !< intermediate def grad at start of FE inc
    crystallite_partionedFi0,&                                                                      !< intermediate def grad at start of homog inc
    crystallite_F0, &                                                                               !< def grad at start of FE inc
    crystallite_partionedF,  &                                                                      !< def grad to be reached at end of homog inc
    crystallite_partionedF0, &                                                                      !< def grad at start of homog inc
    crystallite_Lp, &                                                                               !< current plastic velocitiy grad (end of converged time step)
    crystallite_Lp0, &                                                                              !< plastic velocitiy grad at start of FE inc
    crystallite_partionedLp0, &                                                                     !< plastic velocity grad at start of homog inc
    crystallite_Li, &                                                                               !< current intermediate velocitiy grad (end of converged time step)
    crystallite_Li0, &                                                                              !< intermediate velocitiy grad at start of FE inc
    crystallite_partionedLi0                                                                        !< intermediate velocity grad at start of homog inc
  real(pReal),                dimension(:,:,:,:,:),    allocatable, private :: &
    crystallite_subS0, &                                                                            !< 2nd Piola-Kirchhoff stress vector at start of crystallite inc
    crystallite_invFp, &                                                                            !< inverse of current plastic def grad (end of converged time step)
    crystallite_subFp0,&                                                                            !< plastic def grad at start of crystallite inc
    crystallite_invFi, &                                                                            !< inverse of current intermediate def grad (end of converged time step)
    crystallite_subFi0,&                                                                            !< intermediate def grad at start of crystallite inc
    crystallite_subF,  &                                                                            !< def grad to be reached at end of crystallite inc
    crystallite_subF0, &                                                                            !< def grad at start of crystallite inc
    crystallite_subLp0,&                                                                            !< plastic velocity grad at start of crystallite inc
    crystallite_subLi0                                                                              !< intermediate velocity grad at start of crystallite inc
  real(pReal),                dimension(:,:,:,:,:,:,:), allocatable, public :: &
    crystallite_dPdF                                                                                !< current individual dPdF per grain (end of converged time step)
  logical,                    dimension(:,:,:),         allocatable, public :: &
    crystallite_requested                                                                           !< used by upper level (homogenization) to request crystallite calculation
  logical,                    dimension(:,:,:),         allocatable, private :: &
    crystallite_converged, &                                                                        !< convergence flag
    crystallite_todo, &                                                                             !< flag to indicate need for further computation
    crystallite_localPlasticity                                                                     !< indicates this grain to have purely local constitutive law
 
  enum, bind(c)
    enumerator :: undefined_ID, &
                  phase_ID, &
                  texture_ID, &
                  volume_ID, &
                  orientation_ID, &
                  grainrotation_ID, &
                  defgrad_ID, &
                  fe_ID, &
                  fp_ID, &
                  fi_ID, &
                  lp_ID, &
                  li_ID, &
                  p_ID, &
                  s_ID, &
                  elasmatrix_ID, &
                  neighboringip_ID, &
                  neighboringelement_ID
  end enum
  integer(kind(undefined_ID)),dimension(:,:),   allocatable,          private :: &
    crystallite_outputID                                                                            !< ID of each post result output
    
  type, private :: tOutput                                                                          !< new requested output (per phase)
    character(len=65536), allocatable, dimension(:) :: &
      label
  end type tOutput
  type(tOutput), allocatable, dimension(:), private :: output_constituent
  
  procedure(), pointer :: integrateState
 
  public :: &
    crystallite_init, &
    crystallite_stress, &
    crystallite_stressTangent, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    crystallite_postResults, &
    crystallite_results
  private :: &
    integrateStress, &
    integrateState, &
    integrateStateFPI, &
    integrateStateEuler, &
    integrateStateAdaptiveEuler, &
    integrateStateRK4, &
    integrateStateRKCK45, &
    stateJump

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init
#ifdef DEBUG
  use debug, only: &
    debug_info, &
    debug_reset, &
    debug_level, &
    debug_crystallite, &
    debug_levelBasic
#endif
  use numerics, only: &
    numerics_integrator, &
    worldrank, &
    usePingPong
  use math, only: &
    math_I3, &
    math_EulerToR, &
    math_inv33
  use mesh, only: &
    theMesh, &
    mesh_element
  use IO, only: &
    IO_stringValue, &
    IO_write_jobFile, &
    IO_error
  use material
  use config, only: &
   config_deallocate, &
   config_crystallite, &
   config_phase, &
   crystallite_name
  use constitutive, only: &
    constitutive_initialFi, &
    constitutive_microstructure                                                                     ! derived (shortcut) quantities of given state
 
  implicit none
 
  integer, parameter :: FILEUNIT=434
  logical, dimension(:,:), allocatable :: devNull
  integer :: &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    o = 0, &                                                                                        !< counter in output loop
    r, &  
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax, &                                                                                         !< maximum number of elements
    myNcomponents, &                                                                                !< number of components at current IP
    mySize
 
  character(len=65536), dimension(:), allocatable :: str
 
  write(6,'(/,a)')   ' <<<+-  crystallite init  -+>>>'
 
  cMax = homogenization_maxNgrains
  iMax = theMesh%elem%nIPs
  eMax = theMesh%nElems
 
  allocate(crystallite_S0(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_partionedS0(3,3,cMax,iMax,eMax),       source=0.0_pReal)
  allocate(crystallite_S(3,3,cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_subS0(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_P(3,3,cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_F0(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_partionedF0(3,3,cMax,iMax,eMax),       source=0.0_pReal)
  allocate(crystallite_partionedF(3,3,cMax,iMax,eMax),        source=0.0_pReal)
  allocate(crystallite_subF0(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_subF(3,3,cMax,iMax,eMax),              source=0.0_pReal)
  allocate(crystallite_Fp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedFp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subFp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Fp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_invFp(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_Fi0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedFi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subFi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Fi(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_invFi(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_Fe(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_Lp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedLp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subLp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Lp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_Li0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedLi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subLi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Li(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_dPdF(3,3,3,3,cMax,iMax,eMax),          source=0.0_pReal)
  allocate(crystallite_dt(cMax,iMax,eMax),                    source=0.0_pReal)
  allocate(crystallite_subdt(cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_subFrac(cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_subStep(cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_orientation(cMax,iMax,eMax))
  allocate(crystallite_orientation0(cMax,iMax,eMax))
  allocate(crystallite_localPlasticity(cMax,iMax,eMax),       source=.true.)
  allocate(crystallite_requested(cMax,iMax,eMax),             source=.false.)
  allocate(crystallite_todo(cMax,iMax,eMax),                  source=.false.)
  allocate(crystallite_converged(cMax,iMax,eMax),             source=.true.)
  allocate(crystallite_output(maxval(crystallite_Noutput), &
                              size(config_crystallite))) ;       crystallite_output = ''
  allocate(crystallite_outputID(maxval(crystallite_Noutput), &
                              size(config_crystallite)),         source=undefined_ID)
  allocate(crystallite_sizePostResults(size(config_crystallite)),source=0)
  allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                                      size(config_crystallite)), source=0)
 
  select case(numerics_integrator)
    case(1)
      integrateState => integrateStateFPI
    case(2)
      integrateState => integrateStateEuler
    case(3)
      integrateState => integrateStateAdaptiveEuler
    case(4)
      integrateState => integrateStateRK4
    case(5)
      integrateState => integrateStateRKCK45
  end select
 
 
 
  do c = 1, size(config_crystallite)
#if defined(__GFORTRAN__)
    str = ['GfortranBug86277']
    str = config_crystallite(c)%getStrings('(output)',defaultVal=str)
    if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
    str = config_crystallite(c)%getStrings('(output)',defaultVal=[character(len=65536)::])
#endif
    do o = 1, size(str)
      crystallite_output(o,c) = str(o)
      outputName: select case(str(o))
       case ('phase') outputName
         crystallite_outputID(o,c) = phase_ID
       case ('texture') outputName
         crystallite_outputID(o,c) = texture_ID
       case ('volume') outputName
         crystallite_outputID(o,c) = volume_ID
       case ('orientation') outputName
         crystallite_outputID(o,c) = orientation_ID
       case ('grainrotation') outputName
         crystallite_outputID(o,c) = grainrotation_ID
       case ('defgrad','f') outputName                                                              ! ToDo: no alias (f only)
         crystallite_outputID(o,c) = defgrad_ID
       case ('fe') outputName
         crystallite_outputID(o,c) = fe_ID
       case ('fp') outputName
         crystallite_outputID(o,c) = fp_ID
       case ('fi') outputName
         crystallite_outputID(o,c) = fi_ID
       case ('lp') outputName
         crystallite_outputID(o,c) = lp_ID
       case ('li') outputName
         crystallite_outputID(o,c) = li_ID
       case ('p','firstpiola','1stpiola') outputName                                                ! ToDo: no alias (p only)
         crystallite_outputID(o,c) = p_ID
       case ('s','tstar','secondpiola','2ndpiola') outputName                                       ! ToDo: no alias (s only)
         crystallite_outputID(o,c) = s_ID
       case ('elasmatrix') outputName
         crystallite_outputID(o,c) = elasmatrix_ID
       case ('neighboringip') outputName                                                            ! ToDo: this is not a result, it is static. Should be written out by mesh
         crystallite_outputID(o,c) = neighboringip_ID
       case ('neighboringelement') outputName                                                       ! ToDo: this is not a result, it is static. Should be written out by mesh
         crystallite_outputID(o,c) = neighboringelement_ID
       case default outputName
         call IO_error(105,ext_msg=trim(str(o))//' (Crystallite)')
      end select outputName
    enddo
  enddo
  
  allocate(output_constituent(size(config_phase)))
  do c = 1, size(config_phase)
#if defined(__GFORTRAN__)
    allocate(output_constituent(c)%label(1)) 
    output_constituent(c)%label(1)= 'GfortranBug86277'
    output_constituent(c)%label  = config_phase(c)%getStrings('(output)',defaultVal=output_constituent(c)%label )
    if (output_constituent(c)%label (1) == 'GfortranBug86277') output_constituent(c)%label  = [character(len=pStringLen)::]
#else
    output_constituent(c)%label  = config_phase(c)%getStrings('(output)',defaultVal=[character(len=pStringLen)::])
#endif
  enddo


  do r = 1,size(config_crystallite)
    do o = 1,crystallite_Noutput(r)
      select case(crystallite_outputID(o,r))
        case(phase_ID,texture_ID,volume_ID)
          mySize = 1
        case(orientation_ID,grainrotation_ID)
          mySize = 4
        case(defgrad_ID,fe_ID,fp_ID,fi_ID,lp_ID,li_ID,p_ID,s_ID)
          mySize = 9
        case(elasmatrix_ID)
          mySize = 36
        case(neighboringip_ID,neighboringelement_ID)
          mySize = theMesh%elem%nIPneighbors
        case default
          mySize = 0
      end select
      crystallite_sizePostResult(o,r) = mySize
      crystallite_sizePostResults(r) = crystallite_sizePostResults(r) + mySize
    enddo
  enddo
 
  crystallite_maxSizePostResults = &
    maxval(crystallite_sizePostResults(microstructure_crystallite),microstructure_active)
 

!--------------------------------------------------------------------------------------------------
! write description file for crystallite output
  if (worldrank == 0) then
    call IO_write_jobFile(FILEUNIT,'outputCrystallite')
 
    do r = 1,size(config_crystallite)
      if (any(microstructure_crystallite(mesh_element(4,:)) == r)) then
        write(FILEUNIT,'(/,a,/)') '['//trim(crystallite_name(r))//']'
        do o = 1,crystallite_Noutput(r)
          write(FILEUNIT,'(a,i4)') trim(crystallite_output(o,r))//char(9),crystallite_sizePostResult(o,r)
        enddo
      endif
    enddo
 
    close(FILEUNIT)
  endif
 
  call config_deallocate('material.config/phase')
  call config_deallocate('material.config/crystallite')

!--------------------------------------------------------------------------------------------------
! initialize
 !$OMP PARALLEL DO PRIVATE(myNcomponents,i,c)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNcomponents = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e), FEsolving_execIP(2,e); do c = 1, myNcomponents
      crystallite_Fp0(1:3,1:3,c,i,e) = math_EulerToR(material_EulerAngles(1:3,c,i,e))               ! plastic def gradient reflects init orientation
      crystallite_Fi0(1:3,1:3,c,i,e) = constitutive_initialFi(c,i,e)
      crystallite_F0(1:3,1:3,c,i,e)  = math_I3
      crystallite_localPlasticity(c,i,e) = phase_localPlasticity(material_phase(c,i,e))
      crystallite_Fe(1:3,1:3,c,i,e)  = math_inv33(matmul(crystallite_Fi0(1:3,1:3,c,i,e), &
                                                         crystallite_Fp0(1:3,1:3,c,i,e)))           ! assuming that euler angles are given in internal strain free configuration
      crystallite_Fp(1:3,1:3,c,i,e)  = crystallite_Fp0(1:3,1:3,c,i,e)
      crystallite_Fi(1:3,1:3,c,i,e)  = crystallite_Fi0(1:3,1:3,c,i,e)
      crystallite_requested(c,i,e) = .true.
    enddo; enddo
  enddo
  !$OMP END PARALLEL DO
 
  if(any(.not. crystallite_localPlasticity) .and. .not. usePingPong) call IO_error(601)             ! exit if nonlocal but no ping-pong ToDo: Why not check earlier? or in nonlocal?
 
  crystallite_partionedFp0 = crystallite_Fp0
  crystallite_partionedFi0 = crystallite_Fi0
  crystallite_partionedF0  = crystallite_F0
  crystallite_partionedF   = crystallite_F0
 
  call crystallite_orientations()
  crystallite_orientation0 = crystallite_orientation                                                ! store initial orientations for calculation of grain rotations
 
  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(mesh_element(3,e))
        call constitutive_microstructure(crystallite_Fe(1:3,1:3,c,i,e), &
                                         crystallite_Fp(1:3,1:3,c,i,e), &
                                         c,i,e)                                                     ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
 
  devNull = crystallite_stress()
  call crystallite_stressTangent

#ifdef DEBUG
  if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
    write(6,'(a42,1x,i10)') '    # of elements:                       ', eMax
    write(6,'(a42,1x,i10)') 'max # of integration points/element:     ', iMax
    write(6,'(a42,1x,i10)') 'max # of constituents/integration point: ', cMax
    write(6,'(a42,1x,i10)') 'max # of neigbours/integration point:    ', theMesh%elem%nIPneighbors
    write(6,'(a42,1x,i10)') '    # of nonlocal constituents:          ',count(.not. crystallite_localPlasticity)
    flush(6)
  endif
 
  call debug_info
  call debug_reset
#endif

end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
function crystallite_stress(dummyArgumentToPreventInternalCompilerErrorWithGCC)
  use prec, only: &
    tol_math_check, &
    dNeq0
  use numerics, only: &
    subStepMinCryst, &
    subStepSizeCryst, &
    stepIncreaseCryst
#ifdef DEBUG
  use debug, only: &
    debug_level, &
    debug_crystallite, &
    debug_levelBasic, &
    debug_levelExtensive, &
    debug_levelSelective, &
    debug_e, &
    debug_i, &
    debug_g
#endif
  use IO, only: &
    IO_warning, &
    IO_error
  use math, only: &
    math_inv33
  use mesh, only: &
    theMesh, &
    mesh_element
  use material, only: &
    homogenization_Ngrains, &
    plasticState, &
    sourceState, &
    phase_Nsources, &
    phaseAt, phasememberAt
 
  implicit none
  logical, dimension(theMesh%elem%nIPs,theMesh%Nelems) :: crystallite_stress
  real(pReal), intent(in), optional :: &
    dummyArgumentToPreventInternalCompilerErrorWithGCC
  real(pReal) :: &
    formerSubStep
  integer :: &
    NiterationCrystallite, &                                                                        ! number of iterations in crystallite loop
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    startIP, endIP, &
    s
 
#ifdef DEBUG
  if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0 &
      .and. FEsolving_execElem(1) <= debug_e &
      .and.                          debug_e <= FEsolving_execElem(2)) then
      write(6,'(/,a,i8,1x,i2,1x,i3)')    '<< CRYST stress >> boundary and initial values at el ip ipc ', &
        debug_e,debug_i, debug_g
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> F  ', &
                                          transpose(crystallite_partionedF(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> F0 ', &
                                          transpose(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Fp0', &
                                          transpose(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Fi0', &
                                          transpose(crystallite_partionedFi0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Lp0', &
                                          transpose(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Li0', &
                                          transpose(crystallite_partionedLi0(1:3,1:3,debug_g,debug_i,debug_e))
  endif
#endif

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
  crystallite_subStep = 0.0_pReal
  !$OMP PARALLEL DO
  elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,homogenization_Ngrains(mesh_element(3,e))
      homogenizationRequestsCalculation: if (crystallite_requested(c,i,e)) then
        plasticState    (phaseAt(c,i,e))%subState0(      :,phasememberAt(c,i,e)) = &
        plasticState    (phaseAt(c,i,e))%partionedState0(:,phasememberAt(c,i,e))

        do s = 1, phase_Nsources(phaseAt(c,i,e))
          sourceState(phaseAt(c,i,e))%p(s)%subState0(      :,phasememberAt(c,i,e)) = &
          sourceState(phaseAt(c,i,e))%p(s)%partionedState0(:,phasememberAt(c,i,e))
        enddo
        crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_partionedFp0(1:3,1:3,c,i,e)
        crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_partionedLp0(1:3,1:3,c,i,e)
        crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_partionedFi0(1:3,1:3,c,i,e)
        crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_partionedLi0(1:3,1:3,c,i,e)
        crystallite_subF0(1:3,1:3,c,i,e)  = crystallite_partionedF0(1:3,1:3,c,i,e)
        crystallite_subS0(1:3,1:3,c,i,e)  = crystallite_partionedS0(1:3,1:3,c,i,e)
        crystallite_subFrac(c,i,e) = 0.0_pReal
        crystallite_subStep(c,i,e) = 1.0_pReal/subStepSizeCryst
        crystallite_todo(c,i,e) = .true.
        crystallite_converged(c,i,e) = .false.                                                      ! pretend failed step of 1/subStepSizeCryst
      endif homogenizationRequestsCalculation
    enddo; enddo
  enddo elementLooping1
  !$OMP END PARALLEL DO

  singleRun: if (FEsolving_execELem(1) == FEsolving_execElem(2) .and. &
      FEsolving_execIP(1,FEsolving_execELem(1))==FEsolving_execIP(2,FEsolving_execELem(1))) then
    startIP = FEsolving_execIP(1,FEsolving_execELem(1))
    endIP   = startIP
  else singleRun
    startIP = 1
    endIP = theMesh%elem%nIPs
  endif singleRun

  NiterationCrystallite = 0
  cutbackLooping: do while (any(crystallite_todo(:,startIP:endIP,FEsolving_execELem(1):FEsolving_execElem(2))))
    NiterationCrystallite = NiterationCrystallite + 1

#ifdef DEBUG
    if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0) &
      write(6,'(a,i6)') '<< CRYST stress >> crystallite iteration ',NiterationCrystallite
#endif
    !$OMP PARALLEL DO PRIVATE(formerSubStep)
    elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        do c = 1,homogenization_Ngrains(mesh_element(3,e))
!--------------------------------------------------------------------------------------------------
!  wind forward
          if (crystallite_converged(c,i,e)) then
            formerSubStep = crystallite_subStep(c,i,e)
            crystallite_subFrac(c,i,e) = crystallite_subFrac(c,i,e) + crystallite_subStep(c,i,e)
            crystallite_subStep(c,i,e) = min(1.0_pReal - crystallite_subFrac(c,i,e), &
                                             stepIncreaseCryst * crystallite_subStep(c,i,e))

            crystallite_todo(c,i,e) = crystallite_subStep(c,i,e) > 0.0_pReal                        ! still time left to integrate on?
            if (crystallite_todo(c,i,e)) then
              crystallite_subF0 (1:3,1:3,c,i,e) = crystallite_subF(1:3,1:3,c,i,e)
              crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_Lp  (1:3,1:3,c,i,e)
              crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_Li  (1:3,1:3,c,i,e)
              crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_Fp  (1:3,1:3,c,i,e)
              crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_Fi  (1:3,1:3,c,i,e)
              crystallite_subS0 (1:3,1:3,c,i,e) = crystallite_S   (1:3,1:3,c,i,e)
              !if abbrevation, make c and p private in omp
              plasticState(    phaseAt(c,i,e))%subState0(:,phasememberAt(c,i,e)) &
                = plasticState(phaseAt(c,i,e))%state(    :,phasememberAt(c,i,e))
              do s = 1, phase_Nsources(phaseAt(c,i,e))
                sourceState(    phaseAt(c,i,e))%p(s)%subState0(:,phasememberAt(c,i,e)) &
                  = sourceState(phaseAt(c,i,e))%p(s)%state(    :,phasememberAt(c,i,e))
              enddo
#ifdef DEBUG
              if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0 &
                  .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                         .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
                write(6,'(a,f12.8,a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> winding forward from ', &
                  crystallite_subFrac(c,i,e)-formerSubStep,' to current crystallite_subfrac ', &
                  crystallite_subFrac(c,i,e),' in crystallite_stress at el ip ipc ',e,i,c
#endif
            endif

!--------------------------------------------------------------------------------------------------
!  cut back (reduced time and restore)
          else
            crystallite_subStep(c,i,e)       = subStepSizeCryst * crystallite_subStep(c,i,e)
            crystallite_Fp   (1:3,1:3,c,i,e) =            crystallite_subFp0(1:3,1:3,c,i,e)
            crystallite_invFp(1:3,1:3,c,i,e) = math_inv33(crystallite_Fp    (1:3,1:3,c,i,e))
            crystallite_Fi   (1:3,1:3,c,i,e) =            crystallite_subFi0(1:3,1:3,c,i,e)
            crystallite_invFi(1:3,1:3,c,i,e) = math_inv33(crystallite_Fi    (1:3,1:3,c,i,e))
            crystallite_S    (1:3,1:3,c,i,e) =            crystallite_S0    (1:3,1:3,c,i,e)
            if (crystallite_subStep(c,i,e) < 1.0_pReal) then                                        ! actual (not initial) cutback
              crystallite_Lp (1:3,1:3,c,i,e) =            crystallite_subLp0(1:3,1:3,c,i,e)
              crystallite_Li (1:3,1:3,c,i,e) =            crystallite_subLi0(1:3,1:3,c,i,e)
            endif
            plasticState    (phaseAt(c,i,e))%state(    :,phasememberAt(c,i,e)) &
              = plasticState(phaseAt(c,i,e))%subState0(:,phasememberAt(c,i,e))
            do s = 1, phase_Nsources(phaseAt(c,i,e))
              sourceState(    phaseAt(c,i,e))%p(s)%state(    :,phasememberAt(c,i,e)) &
                = sourceState(phaseAt(c,i,e))%p(s)%subState0(:,phasememberAt(c,i,e))
            enddo

                                                                                                    ! cant restore dotState here, since not yet calculated in first cutback after initialization
            crystallite_todo(c,i,e) = crystallite_subStep(c,i,e) > subStepMinCryst                  ! still on track or already done (beyond repair)
#ifdef DEBUG
            if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
               .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                      .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0)) then
              if (crystallite_todo(c,i,e)) then
                write(6,'(a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> cutback with new crystallite_subStep: ', &
                                                        crystallite_subStep(c,i,e),' at el ip ipc ',e,i,c
              else
                write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> reached minimum step size at el ip ipc ',e,i,c
              endif
            endif
#endif
          endif

!--------------------------------------------------------------------------------------------------
!  prepare for integration
          if (crystallite_todo(c,i,e)) then
            crystallite_subF(1:3,1:3,c,i,e) = crystallite_subF0(1:3,1:3,c,i,e) &
                                            + crystallite_subStep(c,i,e) * (crystallite_partionedF (1:3,1:3,c,i,e) &
                                                                          - crystallite_partionedF0(1:3,1:3,c,i,e))
            crystallite_Fe(1:3,1:3,c,i,e) = matmul(matmul(crystallite_subF (1:3,1:3,c,i,e), &
                                                          crystallite_invFp(1:3,1:3,c,i,e)), &
                                                   crystallite_invFi(1:3,1:3,c,i,e))
            crystallite_subdt(c,i,e) = crystallite_subStep(c,i,e) * crystallite_dt(c,i,e)
            crystallite_converged(c,i,e) = .false.
          endif

        enddo
      enddo
    enddo elementLooping3
    !$OMP END PARALLEL DO

#ifdef DEBUG
    if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0) then
      write(6,'(/,a,f8.5,a,f8.5,/)') '<< CRYST stress >> ',minval(crystallite_subStep), &
                                           ' ≤ subStep ≤ ',maxval(crystallite_subStep)
      write(6,'(/,a,f8.5,a,f8.5,/)') '<< CRYST stress >> ',minval(crystallite_subFrac), &
                                           ' ≤ subFrac ≤ ',maxval(crystallite_subFrac)
      flush(6)
      if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0) then
        write(6,'(/,a,f8.5,1x,a,1x,f8.5,1x,a)') '<< CRYST stress >> subFrac + subStep = ',&
           crystallite_subFrac(debug_g,debug_i,debug_e),'+',crystallite_subStep(debug_g,debug_i,debug_e),'@selective'
        flush(6)
      endif
    endif
#endif
!--------------------------------------------------------------------------------------------------
!  integrate --- requires fully defined state array (basic + dependent state)
    if (any(crystallite_todo)) call integrateState                                                  ! TODO: unroll into proper elementloop to avoid N^2 for single point evaluation
    where(.not. crystallite_converged .and. crystallite_subStep > subStepMinCryst) &                ! do not try non-converged but fully cutbacked any further
      crystallite_todo = .true.                                                                     ! TODO: again unroll this into proper elementloop to avoid N^2 for single point evaluation


  enddo cutbackLooping

! return whether converged or not
  crystallite_stress = .false.
  elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      crystallite_stress(i,e) = all(crystallite_converged(:,i,e)) 
    enddo
  enddo elementLooping5

#ifdef DEBUG
  elementLooping6: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(mesh_element(3,e))
        if (.not. crystallite_converged(c,i,e)) then
          if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
            write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> no convergence at el ip ipc ', &
              e,i,c
        endif
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0)) then
          write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST stress >> solution at el ip ipc ',e,i,c
          write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CRYST stress >> P / MPa', &
                                           transpose(crystallite_P(1:3,1:3,c,i,e))*1.0e-6_pReal
          write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST stress >> Fp', &
                                           transpose(crystallite_Fp(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST stress >> Fi', &
                                           transpose(crystallite_Fi(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST stress >> Lp', &
                                           transpose(crystallite_Lp(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST stress >> Li', &
                                           transpose(crystallite_Li(1:3,1:3,c,i,e))
          flush(6)
        endif
      enddo
    enddo
  enddo elementLooping6
#endif

end function crystallite_stress


!--------------------------------------------------------------------------------------------------
!> @brief calculate tangent (dPdF)
!--------------------------------------------------------------------------------------------------
subroutine crystallite_stressTangent()
  use prec, only: &
    tol_math_check, &
    dNeq0
  use IO, only: &
    IO_warning, &
    IO_error
  use math, only: &
    math_inv33, &
    math_identity2nd, &
    math_3333to99, &
    math_99to3333, &
    math_I3, &
    math_mul3333xx3333, &
    math_mul33xx33, &
    math_invert2, &
    math_det33
  use mesh, only: &
    mesh_element
  use material, only: &
    homogenization_Ngrains
  use constitutive, only:  &
    constitutive_SandItsTangents, &
    constitutive_LpAndItsTangents, &
    constitutive_LiAndItsTangents

  implicit none
  integer :: &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    o, &
    p

  real(pReal), dimension(3,3)     ::   temp_33_1, devNull,invSubFi0, temp_33_2, temp_33_3, temp_33_4
  real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                       dSdF, &
                                       dSdFi, &
                                       dLidS, &
                                       dLidFi, &
                                       dLpdS, &
                                       dLpdFi, &
                                       dFidS, &
                                       dFpinvdF, &
                                       rhs_3333, &
                                       lhs_3333, &
                                       temp_3333
  real(pReal), dimension(9,9)::        temp_99
  logical :: error

  !$OMP PARALLEL DO PRIVATE(dSdF,dSdFe,dSdFi,dLpdS,dLpdFi,dFpinvdF,dLidS,dLidFi,dFidS,invSubFi0,o,p, &
  !$OMP                     rhs_3333,lhs_3333,temp_99,temp_33_1,temp_33_2,temp_33_3,temp_33_4,temp_3333,error)
  elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(mesh_element(3,e))

        call constitutive_SandItsTangents(devNull,dSdFe,dSdFi, &
                                         crystallite_Fe(1:3,1:3,c,i,e), &
                                         crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                       ! call constitutive law to calculate elastic stress tangent
        call constitutive_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e), &
                                           c,i,e)                                                   ! call constitutive law to calculate Li tangent in lattice configuration

        if (sum(abs(dLidS)) < tol_math_check) then
          dFidS = 0.0_pReal
        else
          invSubFi0 = math_inv33(crystallite_subFi0(1:3,1:3,c,i,e))
          lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
          do o=1,3; do p=1,3
            lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                                  + crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidFi(1:3,1:3,o,p))
            lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                                  + crystallite_invFi(1:3,1:3,c,i,e)*crystallite_invFi(p,o,c,i,e)
            rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                                  - crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidS(1:3,1:3,o,p))
          enddo; enddo
          call math_invert2(temp_99,error,math_3333to99(lhs_3333))
          if (error) then
            call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                            ext_msg='inversion error in analytic tangent calculation')
            dFidS = 0.0_pReal
          else
            dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
          endif
          dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
        endif

        call constitutive_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                     ! call constitutive law to calculate Lp tangent in lattice configuration
        dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
        temp_33_1 = transpose(matmul(crystallite_invFp(1:3,1:3,c,i,e), &
                                            crystallite_invFi(1:3,1:3,c,i,e)))
        temp_33_2 = matmul(           crystallite_subF  (1:3,1:3,c,i,e), &
                                  math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)))
        temp_33_3 = matmul(matmul(crystallite_subF  (1:3,1:3,c,i,e), &
                                                crystallite_invFp (1:3,1:3,c,i,e)), &
                                     math_inv33(crystallite_subFi0(1:3,1:3,c,i,e)))

        do o=1,3; do p=1,3
          rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
          temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), &
                                                 crystallite_invFi(1:3,1:3,c,i,e)) &
                                 + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
        enddo; enddo
        lhs_3333 = crystallite_subdt(c,i,e)*math_mul3333xx3333(dSdFe,temp_3333) &
                 + math_mul3333xx3333(dSdFi,dFidS)

        call math_invert2(temp_99,error,math_identity2nd(9)+math_3333to99(lhs_3333))
        if (error) then
          call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                          ext_msg='inversion error in analytic tangent calculation')
          dSdF = rhs_3333
        else
          dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
        endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
        temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
        do o=1,3; do p=1,3
          dFpinvdF(1:3,1:3,p,o) &
            = -crystallite_subdt(c,i,e) &
            * matmul(math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)), &
                            matmul(temp_3333(1:3,1:3,p,o),crystallite_invFi(1:3,1:3,c,i,e)))
        enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
        temp_33_1 = matmul(crystallite_invFp(1:3,1:3,c,i,e), &
                                  matmul(crystallite_S(1:3,1:3,c,i,e), &
                                                transpose(crystallite_invFp(1:3,1:3,c,i,e))))
        temp_33_2 = matmul(crystallite_S(1:3,1:3,c,i,e), &
                                  transpose(crystallite_invFp(1:3,1:3,c,i,e)))
        temp_33_3 = matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                  crystallite_invFp(1:3,1:3,c,i,e))
        temp_33_4 = matmul(matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                                crystallite_invFp(1:3,1:3,c,i,e)), &
                                               crystallite_S(1:3,1:3,c,i,e))

        crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = 0.0_pReal
        do p=1,3
          crystallite_dPdF(p,1:3,p,1:3,c,i,e) = transpose(temp_33_1)
        enddo
        do o=1,3; do p=1,3
          crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
            matmul(matmul(crystallite_subF(1:3,1:3,c,i,e),dFpinvdF(1:3,1:3,p,o)),temp_33_2) + &
            matmul(matmul(temp_33_3,dSdF(1:3,1:3,p,o)),transpose(crystallite_invFp(1:3,1:3,c,i,e))) + &
            matmul(temp_33_4,transpose(dFpinvdF(1:3,1:3,p,o)))
        enddo; enddo

    enddo; enddo
  enddo elementLooping
  !$OMP END PARALLEL DO

end subroutine crystallite_stressTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations
  use math, only: &
    math_rotationalPart33, &
    math_RtoQ
  use material, only: &
    plasticState, &
    material_phase, &
    homogenization_Ngrains
  use mesh, only: &
    mesh_element
  use lattice, only: &
    lattice_qDisorientation
  use plastic_nonlocal, only: &
    plastic_nonlocal_updateCompatibility
 
  implicit none
  integer &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop
 
  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(mesh_element(3,e))
        call crystallite_orientation(c,i,e)%fromRotationMatrix(transpose(math_rotationalPart33(crystallite_Fe(1:3,1:3,c,i,e))))
  enddo; enddo; enddo
  !$OMP END PARALLEL DO
  
  nonlocalPresent: if (any(plasticState%nonLocal)) then
    !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        if (plasticState(material_phase(1,i,e))%nonLocal) &                                         ! if nonlocal model
          call plastic_nonlocal_updateCompatibility(crystallite_orientation,i,e)
    enddo; enddo
    !$OMP END PARALLEL DO
  endif nonlocalPresent

end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(ipc,ip,el, tensor33)
  use math, only: &
   math_inv33, &
   math_EulerToR
  use material, only: &
   material_EulerAngles                                                                             ! ToDo: Why stored? We also have crystallite_orientation0
 
  implicit none
  real(pReal), dimension(3,3) :: crystallite_push33ToRef
  real(pReal), dimension(3,3), intent(in) :: tensor33
  real(pReal), dimension(3,3)             :: T
  integer, intent(in):: &
    el, &
    ip, &
    ipc
 
  T = matmul(math_EulerToR(material_EulerAngles(1:3,ipc,ip,el)), &
                    transpose(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief return results of particular grain
!--------------------------------------------------------------------------------------------------
function crystallite_postResults(ipc, ip, el)
 use math, only: &
   math_qToEuler, &
   math_qToEulerAxisAngle, &
   math_det33, &
   math_I3, &
   inDeg
 use mesh, only: &
   theMesh, &
   mesh_element, &
   mesh_ipVolume, &
   mesh_ipNeighborhood
 use material, only: &
   plasticState, &
   sourceState, &
   microstructure_crystallite, &
   crystallite_Noutput, &
   material_phase, &
   material_texture, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_homogenizedC, &
   constitutive_postResults
 use rotations, only: &
   rotation

 implicit none
 integer, intent(in):: &
   el, &                         !< element index
   ip, &                         !< integration point index
   ipc                           !< grain index

 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,el))) + &
                        1+plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                          sum(sourceState(material_phase(ipc,ip,el))%p(:)%sizePostResults)) :: &
   crystallite_postResults
 real(pReal) :: &
   detF
 integer :: &
   o, &
   c, &
   crystID, &
   mySize, &
   n
 type(rotation) :: rot

 crystID = microstructure_crystallite(mesh_element(4,el))

 crystallite_postResults = 0.0_pReal
 crystallite_postResults(1) = real(crystallite_sizePostResults(crystID),pReal)                  ! header-like information (length)
 c = 1

 do o = 1,crystallite_Noutput(crystID)
   mySize = 0
   select case(crystallite_outputID(o,crystID))
     case (phase_ID)
       mySize = 1
       crystallite_postResults(c+1) = real(material_phase(ipc,ip,el),pReal)                         ! phaseID of grain
     case (texture_ID)
       mySize = 1
       crystallite_postResults(c+1) = real(material_texture(ipc,ip,el),pReal)                       ! textureID of grain
     case (volume_ID)
       mySize = 1
       detF = math_det33(crystallite_partionedF(1:3,1:3,ipc,ip,el))                                 ! V_current = det(F) * V_reference
       crystallite_postResults(c+1) = detF * mesh_ipVolume(ip,el) &
                                           / real(homogenization_Ngrains(mesh_element(3,el)),pReal) ! grain volume (not fraction but absolute)
     case (orientation_ID)
       mySize = 4
       crystallite_postResults(c+1:c+mySize) = crystallite_orientation(ipc,ip,el)%asQuaternion()

     case (grainrotation_ID)
       rot = crystallite_orientation0(ipc,ip,el)%misorientation(crystallite_orientation(ipc,ip,el))
       mySize = 4
       crystallite_postResults(c+1:c+mySize) = rot%asAxisAnglePair()
       crystallite_postResults(c+4) = inDeg * crystallite_postResults(c+4)                          ! angle in degree

! remark: tensor output is of the form 11,12,13, 21,22,23, 31,32,33
! thus row index i is slow, while column index j is fast. reminder: "row is slow"

     case (defgrad_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_partionedF(1:3,1:3,ipc,ip,el)),[mySize])
     case (fe_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fe(1:3,1:3,ipc,ip,el)),[mySize])
     case (fp_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fp(1:3,1:3,ipc,ip,el)),[mySize])
     case (fi_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fi(1:3,1:3,ipc,ip,el)),[mySize])
     case (lp_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Lp(1:3,1:3,ipc,ip,el)),[mySize])
     case (li_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Li(1:3,1:3,ipc,ip,el)),[mySize])
     case (p_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_P(1:3,1:3,ipc,ip,el)),[mySize])
     case (s_ID)
       mySize = 9
       crystallite_postResults(c+1:c+mySize) = &
         reshape(crystallite_S(1:3,1:3,ipc,ip,el),[mySize])
     case (elasmatrix_ID)
       mySize = 36
       crystallite_postResults(c+1:c+mySize) = reshape(constitutive_homogenizedC(ipc,ip,el),[mySize])
     case(neighboringelement_ID)
       mySize = theMesh%elem%nIPneighbors
       crystallite_postResults(c+1:c+mySize) = 0.0_pReal
       forall (n = 1:mySize) &
         crystallite_postResults(c+n) = real(mesh_ipNeighborhood(1,n,ip,el),pReal)
     case(neighboringip_ID)
       mySize = theMesh%elem%nIPneighbors
       crystallite_postResults(c+1:c+mySize) = 0.0_pReal
       forall (n = 1:mySize) &
         crystallite_postResults(c+n) = real(mesh_ipNeighborhood(2,n,ip,el),pReal)
   end select
   c = c + mySize
 enddo

 crystallite_postResults(c+1) = real(plasticState(material_phase(ipc,ip,el))%sizePostResults,pReal)             ! size of constitutive results
 c = c + 1
 if (size(crystallite_postResults)-c > 0) &
   crystallite_postResults(c+1:size(crystallite_postResults)) = &
      constitutive_postResults(crystallite_S(1:3,1:3,ipc,ip,el), crystallite_Fi(1:3,1:3,ipc,ip,el), &
                               ipc, ip, el)

end function crystallite_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine crystallite_results
#if defined(PETSc) || defined(DAMASK_HDF5)
  use lattice
  use results
  use HDF5_utilities
  use rotations
  use config, only: &
    config_name_phase => phase_name                                                                  ! anticipate logical name
   
  use material, only: &
    material_phase_plasticity_type => phase_plasticity
        
  implicit none
  integer :: p,o
  real(pReal),    allocatable, dimension(:,:,:) :: selected_tensors
  type(rotation), allocatable, dimension(:)     :: selected_rotations
  real(pReal),    allocatable, dimension(:,:) :: selected_misorientations 
  character(len=256) :: group,lattice_label
  
  call HDF5_closeGroup(results_addGroup('current/constituent'))   
                                             
  do p=1,size(config_name_phase)
    group = trim('current/constituent')//'/'//trim(config_name_phase(p))
    call HDF5_closeGroup(results_addGroup(group))
    do o = 1, size(output_constituent(p)%label)
      select case (output_constituent(p)%label(o))
        case('f')
          selected_tensors = select_tensors(crystallite_partionedF,p)
          call results_writeDataset(group,selected_tensors,'F',&
                                   'deformation gradient','1')
        case('fe')
          selected_tensors = select_tensors(crystallite_Fe,p)
          call results_writeDataset(group,selected_tensors,'Fe',&
                                   'elastic deformation gradient','1')
        case('fp')
          selected_tensors = select_tensors(crystallite_Fp,p)
          call results_writeDataset(group,selected_tensors,'Fp',&
                                   'plastic deformation gradient','1')
        case('fi')
          selected_tensors = select_tensors(crystallite_Fi,p)
          call results_writeDataset(group,selected_tensors,'Fi',&
                                   'inelastic deformation gradient','1')
        case('lp')
          selected_tensors = select_tensors(crystallite_Lp,p)
          call results_writeDataset(group,selected_tensors,'Lp',&
                                   'plastic velocity gradient','1/s')
        case('li')
          selected_tensors = select_tensors(crystallite_Li,p)
          call results_writeDataset(group,selected_tensors,'Li',&
                                   'inelastic velocity gradient','1/s')
        case('p')
          selected_tensors = select_tensors(crystallite_P,p)
          call results_writeDataset(group,selected_tensors,'P',&
                                   '1st Piola-Kirchoff stress','Pa')
        case('s')
          selected_tensors = select_tensors(crystallite_S,p)
          call results_writeDataset(group,selected_tensors,'S',&
                                   '2nd Piola-Kirchoff stress','Pa')
        case('orientation')
          select case(lattice_structure(p))
            case(LATTICE_iso_ID)
              lattice_label = 'iso'
            case(LATTICE_fcc_ID)
              lattice_label = 'fcc'
            case(LATTICE_bcc_ID)
              lattice_label = 'bcc'
            case(LATTICE_bct_ID)
              lattice_label = 'bct'
            case(LATTICE_hex_ID)
              lattice_label = 'hex'
            case(LATTICE_ort_ID)
              lattice_label = 'ort'
          end select
          selected_rotations = select_rotations(crystallite_orientation,p)
          call results_writeDataset(group,selected_rotations,'orientation',&
                                   'crystal orientation as quaternion',lattice_label)
        
        case('grainrotation')
          selected_misorientations = select_misorientations(crystallite_orientation0,crystallite_orientation,p)                
          call results_writeDataset(group,selected_misorientations,'grain_rotation',&
                                   'grain rotation as axis and angle',lattice_label)
      end select
    enddo
 enddo
 
 contains

!--------------------------------------------------------------------------------------------------
!> @brief select tensors for output
!--------------------------------------------------------------------------------------------------
  function select_tensors(dataset,instance)
 
    use material, only: &
    homogenization_maxNgrains, &
    material_phaseAt
 
    integer, intent(in) :: instance
    real(pReal), dimension(:,:,:,:,:), intent(in) :: dataset
    real(pReal), allocatable, dimension(:,:,:) :: select_tensors
    integer :: e,i,c,j
     
    allocate(select_tensors(3,3,count(material_phaseAt==instance)*homogenization_maxNgrains))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, homogenization_maxNgrains                                                           !ToDo: this needs to be changed for varying Ngrains
        do c = 1, size(material_phaseAt,1)
           if (material_phaseAt(c,e) == instance) then
             j = j + 1
             select_tensors(1:3,1:3,j) = dataset(1:3,1:3,c,i,e)
           endif
        enddo
      enddo
    enddo
   
  end function select_tensors
 
 
!--------------------------------------------------------------------------------------------------
!> @brief select rotations for output
!-------------------------------------------------------------------------------------------------- 
  function select_rotations(dataset,instance)
 
    use material, only: &
    homogenization_maxNgrains, &
    material_phaseAt
 
    integer, intent(in) :: instance
    type(rotation), dimension(:,:,:), intent(in) :: dataset
    type(rotation), allocatable, dimension(:) :: select_rotations
    integer :: e,i,c,j
     
    allocate(select_rotations(count(material_phaseAt==instance)*homogenization_maxNgrains))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, homogenization_maxNgrains                                                            !ToDo: this needs to be changed for varying Ngrains
        do c = 1, size(material_phaseAt,1)
           if (material_phaseAt(c,e) == instance) then
             j = j + 1
             select_rotations(j) = dataset(c,i,e)
           endif
        enddo
      enddo
   enddo
   
 end function select_rotations

!--------------------------------------------------------------------------------------------------
!> @brief select misorientations for output
!-------------------------------------------------------------------------------------------------- 
  function select_misorientations(dataset0,dataset1,instance)
 
    use material, only: &
    homogenization_maxNgrains, &
    material_phaseAt
    use math, only: &
      inDeg        
 
    integer, intent(in) :: instance
    type(rotation), dimension(:,:,:), intent(in) :: dataset0,dataset1
    real(pReal), allocatable, dimension(:,:) :: select_misorientations
    type(rotation) :: rot
    integer :: e,i,c,j
     
    allocate(select_misorientations(4,count(material_phaseAt==instance)*homogenization_maxNgrains))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, homogenization_maxNgrains                                                            !ToDo: this needs to be changed for varying Ngrains
        do c = 1, size(material_phaseAt,1)
           if (material_phaseAt(c,e) == instance) then
             j = j + 1
             rot = dataset0(c,i,e)%misorientation(dataset1(c,i,e))
             select_misorientations(1:4,j) = rot%asAxisAnglePair()
             select_misorientations(4,j) = inDeg * select_misorientations(4,j)
           endif
        enddo
      enddo
   enddo
   
 end function select_misorientations
#endif


end subroutine crystallite_results


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
logical function integrateStress(ipc,ip,el,timeFraction)
  use, intrinsic :: &
    IEEE_arithmetic
  use prec, only:         tol_math_check, &
                          dEq0
  use numerics, only:     nStress, &
                          aTol_crystalliteStress, &
                          rTol_crystalliteStress, &
                          iJacoLpresiduum, &
                          subStepSizeLp, &
                          subStepSizeLi
#ifdef DEBUG
  use debug, only:        debug_level, &
                          debug_e, &
                          debug_i, &
                          debug_g, &
                          debug_crystallite, &
                          debug_levelBasic, &
                          debug_levelExtensive, &
                          debug_levelSelective
#endif

  use constitutive, only: constitutive_LpAndItsTangents, &
                          constitutive_LiAndItsTangents, &
                          constitutive_SandItsTangents
  use math, only:         math_mul33xx33, &
                          math_mul3333xx3333, &
                          math_inv33, &
                          math_det33, &
                          math_I3, &
                          math_identity2nd, &
                          math_3333to99, &
                          math_33to9, &
                          math_9to33
 
  implicit none
  integer, intent(in)::         el, &                                                               ! element index
                                      ip, &                                                         ! integration point index
                                      ipc                                                           ! grain index
  real(pReal), optional, intent(in) :: timeFraction                                                 ! fraction of timestep
 
  real(pReal), dimension(3,3)::       Fg_new, &                                                     ! deformation gradient at end of timestep
                                      Fp_new, &                                                     ! plastic deformation gradient at end of timestep
                                      Fe_new, &                                                     ! elastic deformation gradient at end of timestep
                                      invFp_new, &                                                  ! inverse of Fp_new
                                      Fi_new, &                                                     ! gradient of intermediate deformation stages
                                      invFi_new, &
                                      invFp_current, &                                              ! inverse of Fp_current
                                      invFi_current, &                                              ! inverse of Fp_current
                                      Lpguess, &                                                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                                                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &                                            ! plastic velocity gradient resulting from constitutive law
                                      residuumLp, &                                                 ! current residuum of plastic velocity gradient
                                      residuumLp_old, &                                             ! last residuum of plastic velocity gradient
                                      deltaLp, &                                                    ! direction of next guess
                                      Liguess, &                                                    ! current guess for intermediate velocity gradient
                                      Liguess_old, &                                                ! known last good guess for intermediate velocity gradient
                                      Li_constitutive, &                                            ! intermediate velocity gradient resulting from constitutive law
                                      residuumLi, &                                                 ! current residuum of intermediate velocity gradient
                                      residuumLi_old, &                                             ! last residuum of intermediate velocity gradient
                                      deltaLi, &                                                    ! direction of next guess
                                      S, &                                                          ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                      A, &
                                      B, &
                                      Fe, &                                                         ! elastic deformation gradient
                                      temp_33
  real(pReal), dimension(9) ::        work                                                          ! needed for matrix inversion by LAPACK
  integer,    dimension(9) ::         devNull                                                       ! needed for matrix inversion by LAPACK
  real(pReal), dimension(9,9) ::      dRLp_dLp, &                                                   ! partial derivative of residuum (Jacobian for Newton-Raphson scheme)
                                      dRLp_dLp2, &                                                  ! working copy of dRdLp
                                      dRLi_dLi                                                      ! partial derivative of residuumI (Jacobian for Newton-Raphson scheme)
  real(pReal), dimension(3,3,3,3)::   dS_dFe, &                                                     ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dS_dFi, &
                                      dFe_dLp, &                                                    ! partial derivative of elastic deformation gradient
                                      dFe_dLi, &
                                      dFi_dLi, &
                                      dLp_dFi, &
                                      dLi_dFi, &
                                      dLp_dS, &
                                      dLi_dS
  real(pReal)                         detInvFi, &                                                   ! determinant of InvFi
                                      steplengthLp, &
                                      steplengthLi, &
                                      dt, &                                                         ! time increment
                                      aTolLp, &
                                      aTolLi
  integer                             NiterationStressLp, &                                         ! number of stress integrations
                                      NiterationStressLi, &                                         ! number of inner stress integrations
                                      ierr, &                                                       ! error indicator for LAPACK
                                      o, &
                                      p, &
                                      jacoCounterLp, &
                                      jacoCounterLi                                                 ! counters to check for Jacobian update
  external :: &
    dgesv
 
  !* be pessimistic
  integrateStress = .false.
#ifdef DEBUG
  if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
      .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
  write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> at el ip ipc ',el,ip,ipc
#endif

  if (present(timeFraction)) then
    dt = crystallite_subdt(ipc,ip,el) * timeFraction
    Fg_new = crystallite_subF0(1:3,1:3,ipc,ip,el) &
           + (crystallite_subF(1:3,1:3,ipc,ip,el) - crystallite_subF0(1:3,1:3,ipc,ip,el)) * timeFraction
  else
    dt = crystallite_subdt(ipc,ip,el)
    Fg_new = crystallite_subF(1:3,1:3,ipc,ip,el)
  endif


 !* feed local variables
  Lpguess     =   crystallite_Lp(1:3,1:3,ipc,ip,el)                                                 ! ... and take it as first guess
  Liguess     =   crystallite_Li(1:3,1:3,ipc,ip,el)                                                 ! ... and take it as first guess
  Liguess_old =   Liguess
 
  invFp_current = math_inv33(crystallite_subFp0(1:3,1:3,ipc,ip,el))
  failedInversionFp: if (all(dEq0(invFp_current))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
      write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on inversion of current Fp at el ip ipc ',&
                                     el,ip,ipc
    if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0) &
      write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> current Fp ',transpose(crystallite_subFp0(1:3,1:3,ipc,ip,el))
#endif
    return
  endif failedInversionFp
  A = matmul(Fg_new,invFp_current)                                                            ! intermediate tensor needed later to calculate dFe_dLp

  invFi_current = math_inv33(crystallite_subFi0(1:3,1:3,ipc,ip,el))
  failedInversionFi: if (all(dEq0(invFi_current))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
      write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on inversion of current Fi at el ip ipc ',&
                                     el,ip,ipc
    if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0) &
      write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> current Fi ', &
                                             transpose(crystallite_subFi0(1:3,1:3,ipc,ip,el))
#endif
    return
  endif failedInversionFi
 
  !* start Li loop with normal step length
  NiterationStressLi = 0
  jacoCounterLi      = 0
  steplengthLi       = 1.0_pReal
  residuumLi_old     = 0.0_pReal
 
  LiLoop: do
    NiterationStressLi = NiterationStressLi + 1
    LiLoopLimit: if (NiterationStressLi > nStress) then
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
        write(6,'(a,i3,a,i8,1x,i2,1x,i3,/)') '<< CRYST integrateStress >> reached Li loop limit',nStress, &
                                            ' at el ip ipc ', el,ip,ipc
#endif
      return
    endif LiLoopLimit
 
    invFi_new = matmul(invFi_current,math_I3 - dt*Liguess)
    Fi_new    = math_inv33(invFi_new)
    detInvFi  = math_det33(invFi_new)
 
    !* start Lp loop with normal step length
    NiterationStressLp = 0
    jacoCounterLp      = 0
    steplengthLp       = 1.0_pReal
    residuumLp_old     = 0.0_pReal
    Lpguess_old        = Lpguess
 
    LpLoop: do
      NiterationStressLp = NiterationStressLp + 1
      LpLoopLimit: if (NiterationStressLp > nStress) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
           write(6,'(a,i3,a,i8,1x,i2,1x,i3,/)') '<< CRYST integrateStress >> reached Lp loop limit',nStress, &
                                                ' at el ip ipc ', el,ip,ipc
#endif
        return
      endif LpLoopLimit
 
      !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law
 
      B  = math_I3 - dt*Lpguess
      Fe = matmul(matmul(A,B), invFi_new)
      call constitutive_SandItsTangents(S, dS_dFe, dS_dFi, &
                                        Fe, Fi_new, ipc, ip, el)                                     ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative in unloaded configuration
 
      !* calculate plastic velocity gradient and its tangent from constitutive law
      call constitutive_LpAndItsTangents(Lp_constitutive, dLp_dS, dLp_dFi, &
                                         S, Fi_new, ipc, ip, el)
 
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
        write(6,'(a,i3,/)')                   '<< CRYST integrateStress >> iteration ', NiterationStressLp
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Lpguess', transpose(Lpguess)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Lp_constitutive', transpose(Lp_constitutive)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Fi', transpose(Fi_new)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Fe', transpose(Fe)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> S', transpose(S)
      endif
#endif

      !* update current residuum and check for convergence of loop
      aTolLp = max(rTol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &           ! absolute tolerance from largest acceptable relative error
                   aTol_crystalliteStress)                                                          ! minimum lower cutoff
      residuumLp = Lpguess - Lp_constitutive
 
      if (any(IEEE_is_NaN(residuumLp))) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
          write(6,'(a,i8,1x,i2,1x,i3,a,i3,a)') '<< CRYST integrateStress >> encountered NaN for Lp-residuum at el ip ipc ', &
                                               el,ip,ipc, &
                                               ' ; iteration ', NiterationStressLp,&
                                               ' >> returning..!'
#endif
        return                                                                                      ! ...me = .false. to inform integrator about problem
      elseif (norm2(residuumLp) < aTolLp) then                                                      ! converged if below absolute tolerance
        exit LpLoop                                                                                 ! ...leave iteration loop
      elseif (     NiterationStressLp == 1 &
              .or. norm2(residuumLp) < norm2(residuumLp_old)) then                                  ! not converged, but improved norm of residuum (always proceed in first iteration)...
        residuumLp_old = residuumLp                                                                 ! ...remember old values and...
        Lpguess_old    = Lpguess
        steplengthLp   = 1.0_pReal                                                                  ! ...proceed with normal step length (calculate new search direction)
      else                                                                                          ! not converged and residuum not improved...
        steplengthLp = subStepSizeLp * steplengthLp                                                 ! ...try with smaller step length in same direction
        Lpguess      = Lpguess_old + steplengthLp * deltaLp
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,1x,f7.4)') '<< CRYST integrateStress >> linear search for Lpguess with step', steplengthLp
        endif
#endif
        cycle LpLoop
      endif


     !* calculate Jacobian for correction term
      if (mod(jacoCounterLp, iJacoLpresiduum) == 0) then
        do o=1,3; do p=1,3
          dFe_dLp(o,1:3,p,1:3) = A(o,p)*transpose(invFi_new)                                        ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        enddo; enddo
        dFe_dLp  = - dt * dFe_dLp
        dRLp_dLp = math_identity2nd(9) &
                 - math_3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dS,dS_dFe),dFe_dLp))
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST integrateStress >> dLp_dS', math_3333to99(dLp_dS)
          write(6,'(a,1x,e20.10)')             '<< CRYST integrateStress >> dLp_dS norm', norm2(math_3333to99(dLp_dS))
          write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST integrateStress >> dRLp_dLp', dRLp_dLp-math_identity2nd(9)
          write(6,'(a,1x,e20.10)')             '<< CRYST integrateStress >> dRLp_dLp norm', norm2(dRLp_dLp-math_identity2nd(9))
        endif
#endif
        dRLp_dLp2 = dRLp_dLp                                                                        ! will be overwritten in first call to LAPACK routine
        work = math_33to9(residuumLp)
        call dgesv(9,1,dRLp_dLp2,9,devNull,work,9,ierr)                                             ! solve dRLp/dLp * delta Lp = -res for delta Lp
        if (ierr /= 0) then
#ifdef DEBUG
          if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
           write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on dR/dLp inversion at el ip ipc ', &
                                         el,ip,ipc
            if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
                .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                       .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
              write(6,*)
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dR_dLp',transpose(dRLp_dLp)
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dFe_dLp',transpose(math_3333to99(dFe_dLp))
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dS_dFe (cnst)',transpose(math_3333to99(dS_dFe))
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dLp_dS (cnst)',transpose(math_3333to99(dLp_dS))
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> A',transpose(A)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> B',transpose(B)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Lp_constitutive',transpose(Lp_constitutive)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Lpguess',transpose(Lpguess)
            endif
          endif
#endif
          return
        endif
        deltaLp = - math_9to33(work)
      endif
      jacoCounterLp = jacoCounterLp + 1

      Lpguess = Lpguess + steplengthLp * deltaLp

    enddo LpLoop

   !* calculate intermediate velocity gradient and its tangent from constitutive law
    call constitutive_LiAndItsTangents(Li_constitutive, dLi_dS, dLi_dFi, &
                                       S, Fi_new, ipc, ip, el)

#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
        write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Li_constitutive', transpose(Li_constitutive)
        write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Liguess', transpose(Liguess)
      endif
#endif

    !* update current residuum and check for convergence of loop
    aTolLi = max(rTol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &             ! absolute tolerance from largest acceptable relative error
                 aTol_crystalliteStress)                                                            ! minimum lower cutoff
    residuumLi = Liguess - Li_constitutive
    if (any(IEEE_is_NaN(residuumLi))) then                                                          ! NaN in residuum...
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
        write(6,'(a,i8,1x,i2,1x,i3,a,i3,a)') '<< CRYST integrateStress >> encountered NaN for Li-residuum at el ip ipc ', &
                                             el,ip,ipc, &
                                             ' ; iteration ', NiterationStressLi,&
                                             ' >> returning..!'
#endif
      return                                                                                        ! ...me = .false. to inform integrator about problem
    elseif (norm2(residuumLi) < aTolLi) then                                                        ! converged if below absolute tolerance
      exit LiLoop                                                                                   ! ...leave iteration loop
    elseif (     NiterationStressLi == 1 &
            .or. norm2(residuumLi) < norm2(residuumLi_old)) then                                    ! not converged, but improved norm of residuum (always proceed in first iteration)...
      residuumLi_old = residuumLi                                                                   ! ...remember old values and...
      Liguess_old    = Liguess
      steplengthLi   = 1.0_pReal                                                                    ! ...proceed with normal step length (calculate new search direction)
    else                                                                                            ! not converged and residuum not improved...
      steplengthLi   = subStepSizeLi * steplengthLi                                                 ! ...try with smaller step length in same direction
      Liguess        = Liguess_old + steplengthLi * deltaLi
      cycle LiLoop
    endif
 
    !* calculate Jacobian for correction term
    if (mod(jacoCounterLi, iJacoLpresiduum) == 0) then
      temp_33     = matmul(matmul(A,B),invFi_current)
      do o=1,3; do p=1,3
        dFe_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*temp_33                                             ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        dFi_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*invFi_current
      enddo; enddo
      do o=1,3; do p=1,3
        dFi_dLi(1:3,1:3,o,p) = matmul(matmul(Fi_new,dFi_dLi(1:3,1:3,o,p)),Fi_new)
      enddo; enddo
      dRLi_dLi  = math_identity2nd(9) &
                - math_3333to99(math_mul3333xx3333(dLi_dS, math_mul3333xx3333(dS_dFe, dFe_dLi) + &
                                                  math_mul3333xx3333(dS_dFi, dFi_dLi)))  &
                - math_3333to99(math_mul3333xx3333(dLi_dFi, dFi_dLi))
      work = math_33to9(residuumLi)
      call dgesv(9,1,dRLi_dLi,9,devNull,work,9,ierr)                                                ! solve dRLi/dLp * delta Li = -res for delta Li
      if (ierr /= 0) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
          write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on dR/dLi inversion at el ip ipc ', &
                                        el,ip,ipc
          if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
              .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
            write(6,*)
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dR_dLi',transpose(dRLi_dLi)
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dFe_dLi',transpose(math_3333to99(dFe_dLi))
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dS_dFi (cnst)',transpose(math_3333to99(dS_dFi))
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dLi_dS (cnst)',transpose(math_3333to99(dLi_dS))
            write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Li_constitutive',transpose(Li_constitutive)
            write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Liguess',transpose(Liguess)
          endif
        endif
#endif
        return
      endif
 
      deltaLi = - math_9to33(work)
    endif
    jacoCounterLi = jacoCounterLi + 1
 
    Liguess = Liguess + steplengthLi * deltaLi
  enddo LiLoop
 
  !* calculate new plastic and elastic deformation gradient
  invFp_new = matmul(invFp_current,B)
  invFp_new = invFp_new / math_det33(invFp_new)**(1.0_pReal/3.0_pReal)                              ! regularize
  Fp_new = math_inv33(invFp_new)
  failedInversionInvFp: if (all(dEq0(Fp_new))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on invFp_new inversion at el ip ipc ', &
                                   el,ip,ipc
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
        write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> invFp_new',transpose(invFp_new)
    endif
#endif
    return
  endif failedInversionInvFp
  Fe_new = matmul(matmul(Fg_new,invFp_new),invFi_new)

!--------------------------------------------------------------------------------------------------
! stress integration was successful
  integrateStress = .true.
  crystallite_P    (1:3,1:3,ipc,ip,el) = matmul(matmul(Fg_new,invFp_new),matmul(S,transpose(invFp_new)))
  crystallite_S    (1:3,1:3,ipc,ip,el) = S
  crystallite_Lp   (1:3,1:3,ipc,ip,el) = Lpguess
  crystallite_Li   (1:3,1:3,ipc,ip,el) = Liguess
  crystallite_Fp   (1:3,1:3,ipc,ip,el) = Fp_new
  crystallite_Fi   (1:3,1:3,ipc,ip,el) = Fi_new
  crystallite_Fe   (1:3,1:3,ipc,ip,el) = Fe_new
  crystallite_invFp(1:3,1:3,ipc,ip,el) = invFp_new
  crystallite_invFi(1:3,1:3,ipc,ip,el) = invFi_new

#ifdef DEBUG
  if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0 &
      .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
    write(6,'(a,/)')                     '<< CRYST integrateStress >> successful integration'
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> P / MPa', &
               transpose(crystallite_P(1:3,1:3,ipc,ip,el))*1.0e-6_pReal
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Cauchy / MPa', &
               matmul(crystallite_P(1:3,1:3,ipc,ip,el), transpose(Fg_new)) * 1.0e-6_pReal / math_det33(Fg_new)
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fe Lp Fe^-1', &
               transpose(matmul(Fe_new, matmul(crystallite_Lp(1:3,1:3,ipc,ip,el), math_inv33(Fe_new))))
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fp',transpose(crystallite_Fp(1:3,1:3,ipc,ip,el))
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fi',transpose(crystallite_Fi(1:3,1:3,ipc,ip,el))
  endif
#endif

end function integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
subroutine integrateStateFPI()
#ifdef DEBUG
 use debug, only:        debug_level, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective
#endif
 use numerics, only: &
   nState
 use mesh, only: &
   mesh_element
 use material, only: &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState

 implicit none

 integer :: &
   NiterationState, &                                                                               !< number of iterations in state loop
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s, &
   sizeDotState
 real(pReal) :: &
   zeta
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState) :: &
   residuum_plastic                                                                                 ! residuum for plastic state
 real(pReal), dimension(constitutive_source_maxSizeDotState) :: &
   residuum_source                                                                                  ! residuum for source state
 logical :: &
   doneWithIntegration

 ! --+>> PREGUESS FOR STATE <<+--
 call update_dotState(1.0_pReal)
 call update_state(1.0_pReal)

 NiterationState = 0
 doneWithIntegration = .false.
 crystalliteLooping: do while (.not. doneWithIntegration .and. NiterationState < nState)
   NiterationState = NiterationState + 1

#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0) &
     write(6,'(a,i6)') '<< CRYST stateFPI >> state iteration ',NiterationState
#endif

   ! store previousDotState and previousDotState2
   
   !$OMP PARALLEL DO PRIVATE(p,c)
     do e = FEsolving_execElem(1),FEsolving_execElem(2)
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         do g = 1,homogenization_Ngrains(mesh_element(3,e))
           if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
             p = phaseAt(g,i,e); c = phasememberAt(g,i,e)

             plasticState(p)%previousDotState2(:,c) = merge(plasticState(p)%previousDotState(:,c),&
                                                            0.0_pReal,&
                                                            NiterationState > 1)
             plasticState(p)%previousDotState (:,c) = plasticState(p)%dotState(:,c)
             do s = 1, phase_Nsources(p)
               sourceState(p)%p(s)%previousDotState2(:,c) = merge(sourceState(p)%p(s)%previousDotState(:,c),&
                                                                  0.0_pReal, &
                                                                  NiterationState > 1)
               sourceState(p)%p(s)%previousDotState (:,c) = sourceState(p)%p(s)%dotState(:,c)
             enddo
           endif
       enddo
     enddo
   enddo
   !$OMP END PARALLEL DO

   call update_dependentState
   call update_stress(1.0_pReal)
   call update_dotState(1.0_pReal)
   
   !$OMP PARALLEL
   !$OMP DO PRIVATE(sizeDotState,residuum_plastic,residuum_source,zeta,p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
           p = phaseAt(g,i,e); c = phasememberAt(g,i,e)
           sizeDotState = plasticState(p)%sizeDotState

           zeta = damper(plasticState(p)%dotState         (:,c), &
                         plasticState(p)%previousDotState (:,c), &
                         plasticState(p)%previousDotState2(:,c))
          
           residuum_plastic(1:SizeDotState) = plasticState(p)%state    (1:sizeDotState,c) &
                                            - plasticState(p)%subState0(1:sizeDotState,c)  &
                                            - (  plasticState(p)%dotState        (:,c) * zeta &
                                               + plasticState(p)%previousDotState(:,c) * (1.0_pReal-zeta) &
                                              ) * crystallite_subdt(g,i,e)

           plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%state(1:sizeDotState,c) &
                                                   - residuum_plastic(1:sizeDotState) 
           plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) * zeta &
                                         + plasticState(p)%previousDotState(:,c) * (1.0_pReal - zeta)
           
           crystallite_converged(g,i,e) = converged(residuum_plastic(1:sizeDotState), &
                                                    plasticState(p)%state(1:sizeDotState,c), &
                                                    plasticState(p)%aTolState(1:sizeDotState))
                             

           do s = 1, phase_Nsources(p)
             sizeDotState  = sourceState(p)%p(s)%sizeDotState
             
             zeta = damper(sourceState(p)%p(s)%dotState         (:,c), &
                           sourceState(p)%p(s)%previousDotState (:,c), &
                           sourceState(p)%p(s)%previousDotState2(:,c))

             residuum_source(1:sizeDotState) = sourceState(p)%p(s)%state    (1:sizeDotState,c)  &
                                             - sourceState(p)%p(s)%subState0(1:sizeDotState,c)  &
                                             - (  sourceState(p)%p(s)%dotState         (:,c) * zeta &
                                                 + sourceState(p)%p(s)%previousDotState(:,c) * (1.0_pReal - zeta) &
                                               ) * crystallite_subdt(g,i,e)

             sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%state(1:sizeDotState,c) &
                                                         - residuum_source(1:sizeDotState)
             sourceState(p)%p(s)%dotState(:,c) = sourceState(p)%p(s)%dotState(:,c) * zeta &
                                               + sourceState(p)%p(s)%previousDotState(:,c)* (1.0_pReal - zeta)

             crystallite_converged(g,i,e) = &
             crystallite_converged(g,i,e) .and. converged(residuum_source(1:sizeDotState), &
                                                          sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                                          sourceState(p)%p(s)%aTolState(1:sizeDotState))
           enddo
         endif
   enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. crystallite_converged(g,i,e)) then                                  ! converged and still alive...
         crystallite_todo(g,i,e) = stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP END PARALLEL


   if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck


   ! --- CHECK IF DONE WITH INTEGRATION ---
   doneWithIntegration = .true.
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         doneWithIntegration = .false.
         exit
       endif
     enddo; enddo
   enddo

 enddo crystalliteLooping


 contains

 !--------------------------------------------------------------------------------------------------
 !> @brief calculate the damping for correction of state and dot state
 !--------------------------------------------------------------------------------------------------
 real(pReal) pure function damper(current,previous,previous2)
 
 implicit none
 real(pReal), dimension(:), intent(in) ::&
   current, previous, previous2
 
 real(pReal) :: dot_prod12, dot_prod22
   
 dot_prod12 = dot_product(current  - previous,  previous - previous2)
 dot_prod22 = dot_product(previous - previous2, previous - previous2)
 if ((dot_product(current,previous) < 0.0_pReal .or. dot_prod12 < 0.0_pReal) .and. dot_prod22 > 0.0_pReal) then
   damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
 else
   damper = 1.0_pReal
 endif
   
 end function damper

end subroutine integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief integrate state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
subroutine integrateStateEuler()
 use material, only: &
   plasticState

 implicit none

 call update_dotState(1.0_pReal)
 call update_state(1.0_pReal)
 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck

end subroutine integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
subroutine integrateStateAdaptiveEuler()
 use mesh, only: &
   theMesh, &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState

 implicit none
 integer :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   p, &
   c, &
   s, &
   sizeDotState
   
  ! ToDo: MD: once all constitutives use allocate state, attach residuum arrays to the state in case of adaptive Euler
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems) :: &
   residuum_plastic
 real(pReal), dimension(constitutive_source_maxSizeDotState,&
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems) :: &
   residuum_source

!--------------------------------------------------------------------------------------------------
! contribution to state and relative residui and from Euler integration
 call update_dotState(1.0_pReal)

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,c)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e); c = phasememberAt(g,i,e)
         sizeDotState = plasticState(p)%sizeDotState
         
         residuum_plastic(1:sizeDotState,g,i,e) = plasticState(p)%dotstate(1:sizeDotState,c) &
                                                * (- 0.5_pReal * crystallite_subdt(g,i,e))
         plasticState(p)%state(1:sizeDotState,c) = &
         plasticState(p)%state(1:sizeDotState,c) + plasticState(p)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e) !ToDo: state, partitioned state?
         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
           
           residuum_source(1:sizeDotState,s,g,i,e) = sourceState(p)%p(s)%dotstate(1:sizeDotState,c) &
                                                   * (- 0.5_pReal * crystallite_subdt(g,i,e))
           sourceState(p)%p(s)%state(1:sizeDotState,c) = &
           sourceState(p)%p(s)%state(1:sizeDotState,c) + sourceState(p)%p(s)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e) !ToDo: state, partitioned state?
         enddo
       endif
     enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call update_dotState(1.0_pReal)

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,c)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e); c = phasememberAt(g,i,e)
         sizeDotState = plasticState(p)%sizeDotState
         
         residuum_plastic(1:sizeDotState,g,i,e) = residuum_plastic(1:sizeDotState,g,i,e) &
                                                + 0.5_pReal * plasticState(p)%dotState(:,c) * crystallite_subdt(g,i,e)
              
         crystallite_converged(g,i,e) = converged(residuum_plastic(1:sizeDotState,g,i,e), &
                                                  plasticState(p)%state(1:sizeDotState,c), &
                                                  plasticState(p)%aTolState(1:sizeDotState))

         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
           
           residuum_source(1:sizeDotState,s,g,i,e) = &
           residuum_source(1:sizeDotState,s,g,i,e) + 0.5_pReal * sourceState(p)%p(s)%dotState(:,c) * crystallite_subdt(g,i,e)

           crystallite_converged(g,i,e) = &
           crystallite_converged(g,i,e) .and. converged(residuum_source(1:sizeDotState,s,g,i,e), &
                                                        sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                                        sourceState(p)%p(s)%aTolState(1:sizeDotState))
          enddo
          
       endif
 enddo; enddo; enddo
 !$OMP END PARALLEL DO

 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck
 
end subroutine integrateStateAdaptiveEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 4th order explicit Runge Kutta method
! ToDo: This is totally BROKEN: RK4dotState is never used!!!
!--------------------------------------------------------------------------------------------------
subroutine integrateStateRK4()
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt

 implicit none
 real(pReal), dimension(4), parameter :: &
   TIMESTEPFRACTION = [0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal]                                   ! factor giving the fraction of the original timestep used for Runge Kutta Integration
 real(pReal), dimension(4), parameter :: &
   WEIGHT = [1.0_pReal, 2.0_pReal, 2.0_pReal, 1.0_pReal/6.0_pReal]                                   ! weight of slope used for Runge Kutta integration (final weight divided by 6)

 integer ::                                    e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g, &                                                  ! grain index in grain loop
                                               p, &                                                  ! phase loop
                                               c, &
                                               n, &
                                               s

 call update_dotState(1.0_pReal)


 do n = 1,4

   !$OMP PARALLEL DO PRIVATE(p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e); c = phasememberAt(g,i,e)

         plasticState(p)%RK4dotState(:,c) = WEIGHT(n)*plasticState(p)%dotState(:,c) &
                                          + merge(plasticState(p)%RK4dotState(:,c),0.0_pReal,n>1)
         do s = 1, phase_Nsources(p)
           sourceState(p)%p(s)%RK4dotState(:,c) = WEIGHT(n)*sourceState(p)%p(s)%dotState(:,c) &
                                                + merge(sourceState(p)%p(s)%RK4dotState(:,c),0.0_pReal,n>1)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO

   call update_state(TIMESTEPFRACTION(n))
   call update_deltaState
   call update_dependentState
   call update_stress(TIMESTEPFRACTION(n))
   ! --- dot state and RK dot state---

   first3steps: if (n < 4) then
     call update_dotState(TIMESTEPFRACTION(n))
   endif first3steps

 enddo

 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck

end subroutine integrateStateRK4


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 5th order Runge-Kutta Cash-Karp method with
!> adaptive step size  (use 5th order solution to advance = "local extrapolation")
!--------------------------------------------------------------------------------------------------
subroutine integrateStateRKCK45()
 use mesh, only: &
   mesh_element, &
   theMesh
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState

 implicit none
 real(pReal), dimension(5,5), parameter :: &
   A = reshape([&
     .2_pReal, .075_pReal,   .3_pReal, -11.0_pReal/54.0_pReal,  1631.0_pReal/55296.0_pReal, &
     .0_pReal, .225_pReal,  -.9_pReal,   2.5_pReal,              175.0_pReal/512.0_pReal, &
     .0_pReal,   .0_pReal,  1.2_pReal, -70.0_pReal/27.0_pReal,   575.0_pReal/13824.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,  35.0_pReal/27.0_pReal, 44275.0_pReal/110592.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,    .0_pReal,              253.0_pReal/4096.0_pReal], &
     [5,5], order=[2,1])                                                                            !< coefficients in Butcher tableau (used for preliminary integration in stages 2 to 6)

 real(pReal), dimension(6), parameter :: &
   B = &
     [37.0_pReal/378.0_pReal, .0_pReal, 250.0_pReal/621.0_pReal, &
     125.0_pReal/594.0_pReal, .0_pReal, 512.0_pReal/1771.0_pReal], &                                !< coefficients in Butcher tableau (used for final integration and error estimate)
   DB = B - &
     [2825.0_pReal/27648.0_pReal,    .0_pReal,                18575.0_pReal/48384.0_pReal,&
     13525.0_pReal/55296.0_pReal, 277.0_pReal/14336.0_pReal,      0.25_pReal]                       !< coefficients in Butcher tableau (used for final integration and error estimate)

 real(pReal), dimension(5), parameter :: &
   C = [0.2_pReal, 0.3_pReal, 0.6_pReal, 1.0_pReal, 0.875_pReal]                                    !< coefficients in Butcher tableau (fractions of original time step in stages 2 to 6)

 integer :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   stage, &                                                                                         ! stage index in integration stage loop
   n, &
   p, &
   cc, &
   s, &
   sizeDotState

   ! ToDo: MD: once all constitutives use allocate state, attach residuum arrays to the state in case of RKCK45

 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems) :: &
   residuum_plastic                                                                         ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_source_maxSizeDotState, &
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems) :: &
   residuum_source                                                                   ! relative residuum from evolution in microstructure


 call update_dotState(1.0_pReal)

 ! --- SECOND TO SIXTH RUNGE KUTTA STEP ---

 do stage = 1,5

   ! --- state update ---

   !$OMP PARALLEL DO PRIVATE(p,cc)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e); cc = phasememberAt(g,i,e)

         plasticState(p)%RKCK45dotState(stage,:,cc) = plasticState(p)%dotState(:,cc)
         plasticState(p)%dotState(:,cc) = A(1,stage) * plasticState(p)%RKCK45dotState(1,:,cc)

         do s = 1, phase_Nsources(p)
           sourceState(p)%p(s)%RKCK45dotState(stage,:,cc) = sourceState(p)%p(s)%dotState(:,cc)
           sourceState(p)%p(s)%dotState(:,cc) = A(1,stage) * sourceState(p)%p(s)%RKCK45dotState(1,:,cc)
         enddo

         do n = 2, stage
           plasticState(p)%dotState(:,cc) = plasticState(p)%dotState(:,cc) &
                                          + A(n,stage) * plasticState(p)%RKCK45dotState(n,:,cc)
           do s = 1, phase_Nsources(p)
             sourceState(p)%p(s)%dotState(:,cc) = sourceState(p)%p(s)%dotState(:,cc) &
                                                + A(n,stage) * sourceState(p)%p(s)%RKCK45dotState(n,:,cc)
           enddo
         enddo

       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO

    call update_state(1.0_pReal) !MD: 1.0 correct?
    call update_deltaState
    call update_dependentState
    call update_stress(C(stage))
    call update_dotState(C(stage))

 enddo


!--------------------------------------------------------------------------------------------------
! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE ---

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,cc)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
     if (crystallite_todo(g,i,e)) then
       p = phaseAt(g,i,e); cc = phasememberAt(g,i,e)
       
       sizeDotState = plasticState(p)%sizeDotState
              
       plasticState(p)%RKCK45dotState(6,:,cc) = plasticState (p)%dotState(:,cc)
       
       residuum_plastic(1:sizeDotState,g,i,e) = &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:sizeDotState,cc)),DB) &              ! why transpose? Better to transpose constant DB
       * crystallite_subdt(g,i,e)
       
        plasticState(p)%dotState(:,cc) =  &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:sizeDotState,cc)), B)                ! why transpose? Better to transpose constant B
         
       do s = 1, phase_Nsources(p)
         sizeDotState = sourceState(p)%p(s)%sizeDotState
       
         sourceState(p)%p(s)%RKCK45dotState(6,:,cc) = sourceState(p)%p(s)%dotState(:,cc)
         
         residuum_source(1:sizeDotState,s,g,i,e) = &
           matmul(transpose(sourceState(p)%p(s)%RKCK45dotState(1:6,1:sizeDotState,cc)),DB) &
         * crystallite_subdt(g,i,e)

         sourceState(p)%p(s)%dotState(:,cc)  = &
           matmul(transpose(sourceState(p)%p(s)%RKCK45dotState(1:6,1:sizeDotState,cc)),B)
       enddo
       
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_state(1.0_pReal)
 
 ! --- relative residui and state convergence ---

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,cc)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
       if (crystallite_todo(g,i,e)) then
         p  = phaseAt(g,i,e); cc = phasememberAt(g,i,e)
       
         sizeDotState = plasticState(p)%sizeDotState
         
         crystallite_todo(g,i,e) = converged(residuum_plastic(1:sizeDotState,g,i,e), &
                                             plasticState(p)%state(1:sizeDotState,cc), &
                                             plasticState(p)%aTolState(1:sizeDotState))

         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
         
                  crystallite_todo(g,i,e) = &
                  crystallite_todo(g,i,e) .and. converged(residuum_source(1:sizeDotState,s,g,i,e), &
                                                          sourceState(p)%p(s)%state(1:sizeDotState,cc), &
                                                          sourceState(p)%p(s)%aTolState(1:sizeDotState))
         enddo
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck
 
end subroutine integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief sets convergence flag for nonlocal calculations
!> @detail one non-converged nonlocal sets all other nonlocals to non-converged to trigger cut back
!--------------------------------------------------------------------------------------------------
subroutine nonlocalConvergenceCheck()

 implicit none
 
 if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) &                    ! any non-local not yet converged (or broken)...
   where( .not. crystallite_localPlasticity) crystallite_converged = .false.

end subroutine nonlocalConvergenceCheck


!--------------------------------------------------------------------------------------------------
!> @brief Sets convergence flag based on "todo": every point that survived the integration (todo is
! still .true. is considered as converged
!> @details: For explicitEuler, RK4 and RKCK45, adaptive Euler and FPI have their on criteria
!--------------------------------------------------------------------------------------------------
subroutine setConvergenceFlag()
 use mesh, only: &
   mesh_element
 implicit none
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g                                                                                                !< grain index in grain loop
 
 !OMP DO PARALLEL PRIVATE
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
       crystallite_converged(g,i,e) = crystallite_todo(g,i,e) .or. crystallite_converged(g,i,e)     ! if still "to do" then converged per definition
 enddo; enddo; enddo
 !OMP END DO PARALLEL

end subroutine setConvergenceFlag


 !--------------------------------------------------------------------------------------------------
 !> @brief determines whether a point is converged
 !--------------------------------------------------------------------------------------------------
 logical pure function converged(residuum,state,aTol)
  use prec, only: &
    dEq0
  use numerics, only: &
    rTol => rTol_crystalliteState
    
  implicit none
  real(pReal), intent(in), dimension(:) ::&
    residuum, state, aTol

  converged = all(abs(residuum) <= max(aTol, rTol*abs(state)))

 end function converged


!--------------------------------------------------------------------------------------------------
!> @brief Standard forwarding of state as state = state0 + dotState * (delta t)
!--------------------------------------------------------------------------------------------------
subroutine update_stress(timeFraction)
 use mesh, only: &
   mesh_element
 implicit none
 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g

 !$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_todo(g,i,e) = integrateStress(g,i,e,timeFraction)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then                  ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_stress

!--------------------------------------------------------------------------------------------------
!> @brief tbd
!--------------------------------------------------------------------------------------------------
subroutine update_dependentState()
 use mesh, only: &
   mesh_element
 use constitutive, only: &
   constitutive_dependentState => constitutive_microstructure

 implicit none
 integer ::                                    e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g                                                     ! grain index in grain loop

 !$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(mesh_element(3,e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_dependentState(crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Standard forwarding of state as state = state0 + dotState * (delta t)
!--------------------------------------------------------------------------------------------------
subroutine update_state(timeFraction)
 use material, only: &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use mesh, only: &
   mesh_element

 implicit none
 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s, &
   mySize

 !$OMP PARALLEL DO PRIVATE(mySize,p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       p = phaseAt(g,i,e); c = phasememberAt(g,i,e)

       mySize = plasticState(p)%sizeDotState
       plasticState(p)%state(1:mySize,c) = plasticState(p)%subState0(1:mySize,c) &
                                         + plasticState(p)%dotState (1:mySize,c) &
                                         * crystallite_subdt(g,i,e) * timeFraction
       do s = 1, phase_Nsources(p)
         mySize = sourceState(p)%p(s)%sizeDotState
         sourceState(p)%p(s)%state(1:mySize,c) = sourceState(p)%p(s)%subState0(1:mySize,c) &
                                               + sourceState(p)%p(s)%dotState (1:mySize,c) &
                                               * crystallite_subdt(g,i,e) * timeFraction
       enddo
     endif
 enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_state


!--------------------------------------------------------------------------------------------------
!> @brief triggers calculation of all new rates
!> if NaN occurs, crystallite_todo is set to FALSE. Any NaN in a nonlocal propagates to all others
!--------------------------------------------------------------------------------------------------
subroutine update_dotState(timeFraction)
 use, intrinsic :: &
   IEEE_arithmetic
 use material, only: &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources
 use mesh, only: &
   mesh_element
 use constitutive, only: &
   constitutive_collectDotState

 implicit none
 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s 
 logical :: &
   NaN, &
   nonlocalStop
   
   nonlocalStop = .false.

   !$OMP PARALLEL DO PRIVATE (p,c,NaN)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
         !$OMP FLUSH(nonlocalStop)
        if ((crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) .and. .not. nonlocalStop) then
           call constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                             crystallite_Fe, &
                                             crystallite_Fi(1:3,1:3,g,i,e), &
                                             crystallite_Fp, &
                                             crystallite_subdt(g,i,e)*timeFraction, g,i,e)
           p = phaseAt(g,i,e); c = phasememberAt(g,i,e)
           NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
           do s = 1, phase_Nsources(p)
             NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(s)%dotState(:,c)))
           enddo
           if (NaN) then
             crystallite_todo(g,i,e) = .false.                                                      ! this one done (and broken)
             if (.not. crystallite_localPlasticity(g,i,e)) nonlocalStop = .True.
           endif
         endif
   enddo; enddo; enddo
   !$OMP END PARALLEL DO

 if (nonlocalStop) crystallite_todo = crystallite_todo .and. crystallite_localPlasticity 

end subroutine update_DotState


subroutine update_deltaState
 use, intrinsic :: &
   IEEE_arithmetic
 use prec, only: &
   dNeq0
 use mesh, only: &
   mesh_element
 use material, only: &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use constitutive, only: &
   constitutive_collectDeltaState
 implicit none
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   mySize, &
   myOffset, &
   c, &
   s 
 logical :: &
   NaN, &
   nonlocalStop
   
   nonlocalStop = .false.

   !$OMP PARALLEL DO PRIVATE(p,c,myOffset,mySize,NaN)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
         !$OMP FLUSH(nonlocalStop)
         if ((crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) .and. .not. nonlocalStop) then
        call constitutive_collectDeltaState(crystallite_S(1:3,1:3,g,i,e), &
                                            crystallite_Fe(1:3,1:3,g,i,e), &
                                            crystallite_Fi(1:3,1:3,g,i,e), &
                                            g,i,e)
         p = phaseAt(g,i,e); c = phasememberAt(g,i,e)
         myOffset = plasticState(p)%offsetDeltaState
         mySize   = plasticState(p)%sizeDeltaState
         NaN = any(IEEE_is_NaN(plasticState(p)%deltaState(1:mySize,c)))
         
         if (.not. NaN) then
         
           plasticState(p)%state(myOffset + 1: myOffset + mySize,c) = &
           plasticState(p)%state(myOffset + 1: myOffset + mySize,c) + plasticState(p)%deltaState(1:mySize,c)
           do s = 1, phase_Nsources(p)
             myOffset = sourceState(p)%p(s)%offsetDeltaState
             mySize   = sourceState(p)%p(s)%sizeDeltaState
             NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(s)%deltaState(1:mySize,c)))
             
             if (.not. NaN) then
               sourceState(p)%p(s)%state(myOffset + 1:myOffset + mySize,c) = &
               sourceState(p)%p(s)%state(myOffset + 1:myOffset + mySize,c) + sourceState(p)%p(s)%deltaState(1:mySize,c)
             endif
          enddo
        endif
         
         crystallite_todo(g,i,e) = .not. NaN
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) nonlocalStop = .true.
         endif
       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO
 if (nonlocalStop) crystallite_todo = crystallite_todo .and. crystallite_localPlasticity
 
end subroutine update_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates a jump in the state according to the current state and the current stress
!> returns true, if state jump was successfull or not needed. false indicates NaN in delta state
!--------------------------------------------------------------------------------------------------
logical function stateJump(ipc,ip,el)
 use, intrinsic :: &
   IEEE_arithmetic
 use prec, only: &
   dNeq0
#ifdef DEBUG
 use debug, only: &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_level, &
   debug_crystallite, &
   debug_levelExtensive, &
   debug_levelSelective
#endif
 use material, only: &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use mesh, only: &
   mesh_element
 use constitutive, only: &
   constitutive_collectDeltaState

 implicit none
 integer, intent(in):: &
   el, &                       ! element index
   ip, &                       ! integration point index
   ipc                         ! grain index

 integer :: &
   c, &
   p, &
   mySource, &
   myOffset, &
   mySize

 c = phasememberAt(ipc,ip,el)
 p = phaseAt(ipc,ip,el)

 call constitutive_collectDeltaState(crystallite_S(1:3,1:3,ipc,ip,el), &
                                     crystallite_Fe(1:3,1:3,ipc,ip,el), &
                                     crystallite_Fi(1:3,1:3,ipc,ip,el), &
                                     ipc,ip,el)

 myOffset = plasticState(p)%offsetDeltaState
 mySize   = plasticState(p)%sizeDeltaState

 if( any(IEEE_is_NaN(plasticState(p)%deltaState(1:mySize,c)))) then                                       ! NaN occured in deltaState
   stateJump = .false.
   return
 endif

 plasticState(p)%state(myOffset + 1:myOffset + mySize,c) = &
 plasticState(p)%state(myOffset + 1:myOffset + mySize,c) + plasticState(p)%deltaState(1:mySize,c)

 do mySource = 1, phase_Nsources(p)
   myOffset = sourceState(p)%p(mySource)%offsetDeltaState
   mySize   = sourceState(p)%p(mySource)%sizeDeltaState
   if (any(IEEE_is_NaN(sourceState(p)%p(mySource)%deltaState(1:mySize,c)))) then   ! NaN occured in deltaState
     stateJump = .false.
     return
   endif
   sourceState(p)%p(mySource)%state(myOffset + 1: myOffset + mySize,c) = &
   sourceState(p)%p(mySource)%state(myOffset + 1: myOffset + mySize,c) + sourceState(p)%p(mySource)%deltaState(1:mySize,c)
 enddo

#ifdef DEBUG
 if (any(dNeq0(plasticState(p)%deltaState(1:mySize,c))) &
     .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
   write(6,'(a,i8,1x,i2,1x,i3, /)') '<< CRYST >> update state at el ip ipc ',el,ip,ipc
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> deltaState', plasticState(p)%deltaState(1:mySize,c)
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', &
     plasticState(p)%state(myOffset + 1                : &
                           myOffset + mySize,c)
 endif
#endif

 stateJump = .true.

end function stateJump

end module crystallite
