! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_amps_diagnostics

  use mpas_derived_types, only : MPAS_pool_type
  use mpas_derived_types, only : MPAS_clock_type
  use mpas_derived_types, only : MPAS_LOG_ERR
  use mpas_derived_types, only : MPAS_LOG_CRIT
  use mpas_log          , only : mpas_log_write
  use mpas_kind_types   , only : RKIND
  use mpas_constants    , only : RVORD
  use mpas_constants    , only : RGAS

  type (MPAS_pool_type) , pointer :: mesh
  type (MPAS_clock_type), pointer :: clock
  type (MPAS_pool_type) , pointer :: state
  type (MPAS_pool_type) , pointer :: diag
  type (MPAS_clock_type), pointer :: simulationClock

  public :: amps_diagnostics_setup, &
       amps_diagnostics_update, &
       amps_diagnostics_compute, &
       amps_diagnostics_reset, &
       amps_diagnostics_cleanup

  private


contains

  !-----------------------------------------------------------------------
  !  routine amps_diagnostics_setup
  !
  !> \brief Initialize the diagnostic
  !> \author 
  !> \date   
  !> \details
  !>  Initialize the diagnostic
  !
  !-----------------------------------------------------------------------
  subroutine amps_diagnostics_setup(all_pools, simulation_clock)

    use mpas_derived_types, only : MPAS_pool_type
    use mpas_derived_types, only : MPAS_clock_type
    use mpas_pool_routines, only : mpas_pool_get_subpool

    implicit none

    type (MPAS_pool_type) , pointer :: all_pools
    type (MPAS_clock_type), pointer :: simulation_clock

    call mpas_log_write("This begins processing in AMPS_DIAGNOSTICS_SETUP.")

    ! Perform initialization, memory allocation, etc.

    ! Also, save pointers to any pools that will be used by this diagnostic
    ! E.g.,
    call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
    call mpas_pool_get_subpool(all_pools, 'state', state)
    call mpas_pool_get_subpool(all_pools, 'diag', diag)

    clock => simulation_clock
    call mpas_log_write("This completes processing in AMPS_DIAGNOSTICS_SETUP.")

  end subroutine amps_diagnostics_setup


  !-----------------------------------------------------------------------
  !  routine amps_diagnostics_update
  !
  !> \brief Handle diagnostic calculation needed after each timestep
  !> \author 
  !> \date   
  !> \details
  !>  Handle diagnostic calculation needed after each timestep
  !
  !-----------------------------------------------------------------------
  subroutine amps_diagnostics_update()

    implicit none

    ! Called at the end of every timestep
    ! Update extrema, accumulations, etc.

    call mpas_log_write("This a call to AMPS_DIAGNOSTICS_UPDATE, which currently does nothing.")

  end subroutine amps_diagnostics_update


  !-----------------------------------------------------------------------
  !  routine amps_diagnostics_compute
  !
  !> \brief Compute diagnostic before model output is written
  !> \author 
  !> \date   
  !> \details
  !>  Compute diagnostic before model output is written
  !
  !-----------------------------------------------------------------------
  subroutine amps_diagnostics_compute()

    ! Called immediately before diagnostics will be written.
    ! Compute the diagnostic

    use mpas_pool_routines, only : mpas_pool_get_dimension
    use mpas_pool_routines, only : mpas_pool_get_array
    use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

    implicit none

    integer, pointer :: nCells
    integer, pointer :: nVertLevels
    ! Fields that are computed in this routine
    real (kind=RKIND), dimension(:), pointer :: R4ceiling

    ! Other fields used in the computation of amps diagnostics 
    ! defined above
    integer, pointer :: index_qv
    integer, pointer :: index_qc
    integer, pointer :: index_qr
    integer, pointer :: index_qi
    integer, pointer :: index_qs
    real (kind=RKIND), dimension(:,:)  , pointer :: zgrid
    real (kind=RKIND), dimension(:,:)  , pointer :: exner
    real (kind=RKIND), dimension(:,:)  , pointer :: pressure_p
    real (kind=RKIND), dimension(:,:)  , pointer :: pressure_base
    real (kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=RKIND), dimension(:,:)  , pointer :: theta_m
    real (kind=RKIND), allocatable, dimension(:) :: p8w
    real (kind=RKIND), allocatable, dimension(:) :: p

    real (kind=RKIND), parameter :: CLGFAC = 1.0
    real (kind=RKIND), parameter :: THCON = 0.02
    real (kind=RKIND), parameter :: EPS = 0.622
    real (kind=RKIND), parameter :: RHOWAT = 1.E3
    real (kind=RKIND), parameter :: RHOICE = 917.0
    real (kind=RKIND), parameter :: COEFLC = 144.7
    real (kind=RKIND), parameter :: COEFLP = 2.24
    real (kind=RKIND), parameter :: COEFFC = 327.8
    real (kind=RKIND), parameter :: COEFFP = CLGFAC*10.36
    real (kind=RKIND), parameter :: EXPONLC = 0.8800
    real (kind=RKIND), parameter :: EXPONLP = 0.7500
    real (kind=RKIND), parameter :: EXPONFC = 1.0000
    real (kind=RKIND), parameter :: EXPONFP = 0.7776
    
    integer :: iCell
    integer :: k
    real (kind=RKIND) :: exting
    real (kind=RKIND) :: pressure
    real (kind=RKIND) :: temperature
    real (kind=RKIND) :: t_virtual
    real (kind=RKIND) :: rho_air
    real (kind=RKIND) :: vovrmd
    real (kind=RKIND) :: conclc
    real (kind=RKIND) :: conclp
    real (kind=RKIND) :: concfc
    real (kind=RKIND) :: concfp
    real (kind=RKIND) :: beta
    real (kind=RKIND) :: dz
    real (kind=RKIND) :: extingt
    

    if (.not. MPAS_field_will_be_written("ceiling_RIP4")) return

    call mpas_log_write("This begins processing in AMPS_DIAGNOSTICS_COMPUTE.")

    call mpas_pool_get_dimension(mesh, 'nCells', nCells)
    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(state, 'index_qv', index_qv)
    call mpas_pool_get_dimension(state, 'index_qc', index_qc)
    call mpas_pool_get_dimension(state, 'index_qi', index_qi)
    call mpas_pool_get_dimension(state, 'index_qr', index_qr)
    call mpas_pool_get_dimension(state, 'index_qs', index_qs)
    call mpas_pool_get_array(mesh, 'zgrid', zgrid)
    call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
    call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
    call mpas_pool_get_array(diag, 'exner', exner)
    call mpas_pool_get_array(diag, 'ceiling_RIP4', R4ceiling)
    call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
    call mpas_pool_get_array(state, 'scalars', scalars, 1)

    !  Calculate

    CellLoop : do iCell = 1, nCells

        exting = 1.0
        R4ceiling(iCell) = 0.0
        KLOOP : do k = 1, nVertLevels
            pressure = (pressure_base(k,iCell) + pressure_p(k,iCell))
            temperature = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell)
            t_virtual = temperature * (EPS+scalars(index_qv,k,iCell))/(EPS*(1+scalars(index_qv,k,iCell)))
            rho_air = pressure/(RGAS*t_virtual)
            vovrmd = (1+scalars(index_qv,k,iCell))/rho_air &
                 &  + (scalars(index_qc,k,iCell)+scalars(index_qr,k,iCell))/RHOWAT &
                 &  + (scalars(index_qi,k,iCell)+scalars(index_qs,k,iCell))/RHOICE
            conclc = 1.E3*scalars(index_qc,k,iCell)/vovrmd
            conclp = 1.E3*scalars(index_qr,k,iCell)/vovrmd
            concfc = 1.E3*scalars(index_qi,k,iCell)/vovrmd
            concfp = 1.E3*scalars(index_qs,k,iCell)/vovrmd
            beta = COEFFC * concfc**EXPONFC + COEFFP * concfp**EXPONFP &
                 &   +  COEFLC * conclc**EXPONLC + COEFLP * conclp**EXPONLP &
                 &   +  1.e-10
            dz = 1.E-3 * (zgrid(k+1,iCell) - zgrid(k,iCell))  ! 1.E-3 term to convert from meters to km

            !
            !  Integrate up from surface.  extingt is 1.0 at the (full-level) surface,
            !  and is reduced for the microphysical species within a layer as you go up.
            !
            
            extingt = exting * exp(-beta*dz)
            if (extingt <= THCON) then
                R4ceiling(iCell) = R4ceiling(iCell) - log(THCON/exting)/beta
                exit KLOOP
            endif
            R4ceiling(iCell) = R4ceiling(iCell) + dz
            exting = extingt
        enddo KLOOP
        R4ceiling(iCell) = min(R4ceiling(iCell)*1000., real(9000,kind=RKIND))
    enddo CellLoop

  end subroutine amps_diagnostics_compute

  !-----------------------------------------------------------------------
  !  routine amps_diagnostics_reset
  !
  !> \brief Reset diagnostic after it has been written
  !> \author 
  !> \date   
  !> \details
  !>  Reset diagnostic after it has been written
  !
  !-----------------------------------------------------------------------
  subroutine amps_diagnostics_reset()

    implicit none

    ! Called immediately after diagnostics have been written
    ! Reset counters, accumulations, etc.
    call mpas_log_write("This a call to AMPS_DIAGNOSTICS_RESET, which currently does nothing.")

  end subroutine amps_diagnostics_reset


  !-----------------------------------------------------------------------
  !  routine amps_diagnostics_cleanup
  !
  !> \brief Finalizes diagnostic
  !> \author Michael Duda
  !> \date   6 September 2016
  !> \details
  !>  Finalizes diagnostic
  !
  !-----------------------------------------------------------------------
  subroutine amps_diagnostics_cleanup()

    implicit none

    ! Deallocate scratch arrays, etc.
    call mpas_log_write("This a call to AMPS_DIAGNOSTICS_CLEANUP, which currently does nothing.")

  end subroutine amps_diagnostics_cleanup

end module mpas_amps_diagnostics
