! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_crtm_diagnostics

  use mpas_derived_types, only : MPAS_pool_type
  use mpas_derived_types, only : MPAS_clock_type
  use mpas_derived_types, only : MPAS_time_type
  use mpas_derived_types, only : MPAS_LOG_ERR
  use mpas_derived_types, only : MPAS_LOG_CRIT
  use mpas_log          , only : mpas_log_write
  use mpas_kind_types   , only : StrKIND
  use mpas_kind_types   , only : RKIND
  use mpas_kind_types   , only : R8KIND
  use mpas_constants    , only : rvord
  use mpas_constants    , only : pii
  use mpas_atmphys_constants , only : degrad

  type (MPAS_pool_type) , pointer :: mesh
  type (MPAS_clock_type), pointer :: clock
  type (MPAS_pool_type) , pointer :: state
  type (MPAS_pool_type) , pointer :: diag
  type (MPAS_pool_type) , pointer :: diag_physics
  type (MPAS_pool_type) , pointer :: sfc_input
  type (MPAS_clock_type), pointer :: simulationClock

  character(len=StrKIND), pointer :: config_microp_scheme

  public :: crtm_diagnostics_setup,   &
       &    crtm_diagnostics_update,  &
       &    crtm_diagnostics_compute, &
       &    crtm_diagnostics_reset,   &
       &    crtm_diagnostics_cleanup

  private

contains

  !-----------------------------------------------------------------------
  !  routine crtm_diagnostics_setup
  !
  !> \brief Initialize the diagnostic
  !> \author 
  !> \date   
  !> \details
  !>  Initialize the diagnostic
  !
  !-----------------------------------------------------------------------
  subroutine crtm_diagnostics_setup(configs, all_pools, simulation_clock)

    use mpas_derived_types, only : MPAS_pool_type
    use mpas_derived_types, only : MPAS_clock_type
    use mpas_pool_routines, only : mpas_pool_get_config
    use mpas_pool_routines, only : mpas_pool_get_subpool

    implicit none

    type (MPAS_pool_type) , pointer :: configs
    type (MPAS_pool_type) , pointer :: all_pools
    type (MPAS_clock_type), pointer :: simulation_clock

    call mpas_log_write("This begins processing in CRTM_DIAGNOSTICS_SETUP.")

    ! Perform initialization, memory allocation, etc.

    ! Also, save pointers to any pools that will be used by this diagnostic
    ! E.g.,
    call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
    call mpas_pool_get_subpool(all_pools, 'state', state)
    call mpas_pool_get_subpool(all_pools, 'diag', diag)
    call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
    call mpas_pool_get_subpool(all_pools, 'sfc_input', sfc_input)
    call mpas_pool_get_config(configs, 'config_microp_scheme', config_microp_scheme)

    if (trim(config_microp_scheme) == "mp_wsm6") then
    elseif (trim(config_microp_scheme) == "mp_thompson") then
    else
       call mpas_log_write("CRTM_DIAGNOSTICS_SETUP:  Microphysics scheme unrecognized for CRTM diagnostics: "//trim(config_microp_scheme), messageType=MPAS_LOG_CRIT)
    endif

    clock => simulation_clock
    call mpas_log_write("This completes processing in CRTM_DIAGNOSTICS_SETUP.")

  end subroutine crtm_diagnostics_setup


  !-----------------------------------------------------------------------
  !  routine crtm_diagnostics_update
  !
  !> \brief Handle diagnostic calculation needed after each timestep
  !> \author 
  !> \date   
  !> \details
  !>  Handle diagnostic calculation needed after each timestep
  !
  !-----------------------------------------------------------------------
  subroutine crtm_diagnostics_update()

    implicit none

    ! Called at the end of every timestep
    ! Update extrema, accumulations, etc.

    call mpas_log_write("This a call to CRTM_DIAGNOSTICS_UPDATE, which currently does nothing.")

  end subroutine crtm_diagnostics_update
!-----------------------------------------------------------------------
  !  routine crtm_diagnostics_compute
  !
  !> \brief Compute diagnostic before model output is written
  !> \author 
  !> \date   
  !> \details
  !>  Compute diagnostic before model output is written
  !
  !-----------------------------------------------------------------------
  subroutine crtm_diagnostics_compute()

    ! Called immediately before diagnostics will be written.
    ! Compute the diagnostic

    use mpas_pool_routines, only : mpas_pool_get_dimension
    use mpas_pool_routines, only : mpas_pool_get_array
    use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
    use crtm_module, only: SUCCESS
    use crtm_module, only: RT_ADA
    use crtm_module, only: crtm_channelinfo_type
    use crtm_module, only: crtm_geometry_type
    use crtm_module, only: crtm_atmosphere_type
    use crtm_module, only: crtm_surface_type
    use crtm_module, only: crtm_options_type
    use crtm_module, only: crtm_rtsolution_type
    use crtm_module, only: crtm_init
    use crtm_module, only: crtm_channelinfo_subset
    use crtm_module, only: crtm_channelinfo_n_channels
    use crtm_module, only: crtm_atmosphere_create
    use crtm_module, only: crtm_atmosphere_associated
    use crtm_module, only: crtm_rtsolution_associated
    use crtm_module, only: crtm_rtsolution_create
    use crtm_module, only: crtm_options_create
    use crtm_module, only: crtm_options_setvalue
    use crtm_module, only: crtm_forward
    use crtm_module, only: crtm_destroy
    use crtm_module, only: crtm_channelinfo_destroy
    use crtm_module, only: crtm_atmosphere_destroy
    use crtm_module, only: crtm_geometry_destroy
    use crtm_module, only: crtm_surface_destroy
    use crtm_module, only: crtm_rtsolution_destroy
    use crtm_module, only: crtm_options_destroy

    implicit none

    integer, pointer :: nCells
    integer, pointer :: nVertLevels

    ! Fields that are computed in this routine
    real (kind=RKIND), dimension(:), pointer :: diagnostic_array

    ! Other fields used in the computation of crtm diagnostics 
    ! defined above
    integer :: iCell
    integer :: n
    integer :: n_clouds
    integer :: n_absorbers
    integer :: n_aerosols
    integer :: ierr
    integer :: allocation_status
    integer :: n_channels
    integer :: n_sensors   ! grand total across all satellites
    integer :: nVisible
    integer :: i
    integer :: k
    integer :: idx
    integer, allocatable, dimension(:) :: visible_idx
    type(crtm_geometry_type)   , allocatable, dimension(:) :: geo_subset   ! dimensioned (nVisible)
    type(crtm_atmosphere_type) , allocatable, dimension(:) :: atm          ! dimensioned (nCells)
    type(crtm_atmosphere_type) , allocatable, dimension(:) :: atm_subset   ! dimensioned (nVisible)
    type(crtm_surface_type)    , allocatable, dimension(:) :: sfc          ! dimensioned (nCells)
    type(crtm_surface_type)    , allocatable, dimension(:) :: sfc_subset   ! dimensioned (nVisible)
    type(crtm_rtsolution_type) , allocatable, dimension(:,:) :: rts_subset ! dimensioned (n_channels,nVisible)
    type(crtm_options_type)    , allocatable, dimension(:) :: opt          ! dimensioned (nCells)

    character(len=2048) :: COEFFICIENT_PATH
    character(len=2048) :: IRiceCoeff_File
    character(len=2048) :: IRlandCoeff_File
    character(len=2048) :: IRsnowCoeff_File
    character(len=2048) :: IRwaterCoeff_File
    character(len=2048) :: MWwaterCoeff_File
    character(len=2048) :: VISiceCoeff_File
    character(len=2048) :: VISlandCoeff_File
    character(len=2048) :: VISsnowCoeff_File
    character(len=2048) :: VISwaterCoeff_File

    integer, parameter :: known_diagnostic_count  = 31
    character(len=64) :: diagnostic_field_category
    character(len=64) :: diagnostic_field_name
    character(len=64), dimension(known_diagnostic_count), parameter :: known_diagnostic_field = (/ &
         &  "goes18_radiance_channel_01       ",  &
         &  "goes18_radiance_channel_02       ",  &
         &  "goes18_brightnessT_channel_08    ",  &
         &  "goes18_brightnessT_channel_09    ",  &
         &  "goes18_brightnessT_channel_10    ",  &
         &  "goes19_radiance_channel_01       ",  &
         &  "goes19_radiance_channel_02       ",  &
         &  "goes19_brightnessT_channel_08    ",  &
         &  "goes19_brightnessT_channel_09    ",  &
         &  "goes19_brightnessT_channel_10    ",  &
         &  "himawari9_radiance_channel_01    ",  &
         &  "himawari9_radiance_channel_02    ",  &
         &  "himawari9_radiance_channel_03    ",  &
         &  "himawari9_brightnessT_channel_08 ",  &
         &  "himawari9_brightnessT_channel_09 ",  &
         &  "himawari9_brightnessT_channel_10 ",  &
         &  "meteosat9_radiance_channel_01    ",  &
         &  "meteosat9_radiance_channel_02    ",  &
         &  "meteosat9_radiance_channel_12    ",  &
         &  "meteosat9_brightnessT_channel_05 ",  &
         &  "meteosat9_brightnessT_channel_06 ",  &
         &  "meteosat10_radiance_channel_01   ",  &
         &  "meteosat10_radiance_channel_02   ",  &
         &  "meteosat10_radiance_channel_12   ",  &
         &  "meteosat10_brightnessT_channel_05",  &
         &  "meteosat10_brightnessT_channel_06",  &
         &  "meteosat11_radiance_channel_01   ",  &
         &  "meteosat11_radiance_channel_02   ",  &
         &  "meteosat11_radiance_channel_12   ",  &
         &  "meteosat11_brightnessT_channel_05",  &
         &  "meteosat11_brightnessT_channel_06"   /)

    integer, parameter :: max_channels = 25
    logical, dimension(max_channels) :: flag_vis_channel
    logical, dimension(max_channels) :: flag_ir_channel
    character(len=24) :: visible_sensor_designation
    character(len=24) :: ir_sensor_designation

    integer, parameter :: satellite_count = 6
    character(len=24), dimension(satellite_count), parameter :: satellite_id = (/ &
         &   "goes18    ", &
         &   "goes19    ", &
         &   "himawari9 ", &
         &   "meteosat9 ", &
         &   "meteosat10", &
         &   "meteosat11"  /)
    integer :: idiagnostic
    integer :: isatellite
    integer :: ichannel

    !
    !  One satellite's worth of CHANNEL SELECTION -- NOT ChannelInfo
    !  itself. chinfo/sensor_id deliberately live in ONE persistent flat
    !  array each (see below), never copied into or out of this type,
    !  since CRTM_ChannelInfo_type may hold pointer (not allocatable)
    !  internal components -- a derived-type assignment copy of it,
    !  followed by deallocating the source, would leave a dangling
    !  pointer. channel_subset/nchannels_subset are plain integer data
    !  with no such risk, so they're the only things this type holds.
    !
    type :: satellite_sensors_type
       integer :: n_sensors = 0
       ! Index 1 of channel_subset is the index for the sensor for the satellite in question.
       ! Index 2 of channel_subset is the index for the channel of the sensor referred to in Index 1.
       integer, dimension(10,max_channels) :: channel_subset
       integer, dimension(10) :: nchannels_subset
    end type satellite_sensors_type

    type(satellite_sensors_type) :: sat(satellite_count)
    integer, dimension(satellite_count) :: satellite_chinfo_start

    !
    !  sensor_id/chinfo: ONE persistent flat array each, spanning every
    !  sensor across every satellite, allocated once n_sensors is known
    !  and never copied elsewhere. chinfo in particular is populated in
    !  place by CRTM_Init and stays exactly where Init put it for the
    !  rest of the routine's lifetime.
    !
    character(len=20), allocatable :: sensor_id(:)
    type(crtm_channelinfo_type), allocatable :: chinfo(:)

    !
    !  Per-satellite prescan scratch -- sensor_id only (chinfo is not
    !  needed here at all: setup_sensor_channel_subset never populates
    !  chinfo with anything, it's pure output allocation, so there is
    !  nothing gained by producing it during the prescan).
    !
    character(len=20), allocatable :: local_sensor_id(:)
    integer :: local_n_sensors

    !
    !  Early return-check, to see if we need to calculate anything in this routine.
    !

    do idiagnostic = 1, known_diagnostic_count
       if (MPAS_field_will_be_written(known_diagnostic_field(idiagnostic))) then
          exit
       endif
       if (idiagnostic == known_diagnostic_count) then
          return
       endif
    enddo

    call mpas_log_write("This begins processing in CRTM_DIAGNOSTICS_COMPUTE.")

    if (trim(config_microp_scheme) == "mp_wsm6") then
    else if (trim(config_microp_scheme) == "mp_thompson") then
    else
       call mpas_log_write("Microphysics scheme unrecognized for CRTM diagnostics", messageType=MPAS_LOG_CRIT)
    endif

    call mpas_pool_get_dimension(mesh, 'nCells', nCells)
    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

    n_absorbers = 2
    if (trim(config_microp_scheme) == "mp_wsm6") then
       n_clouds = 5
    else if (trim(config_microp_scheme) == "mp_thompson") then
       n_clouds = 5
    endif
    n_aerosols = 0

    !
    !  Scan all satellites once. Only sensor_id (a plain string, safe to
    !  accumulate) and the per-satellite channel_subset/nchannels_subset
    !  come out of this pass -- chinfo doesn't exist yet at all.
    !

    n_sensors = 0
    allocate(sensor_id(2*satellite_count))

    PRESCAN_SATELLITE_LOOP : do isatellite = 1, satellite_count
       flag_vis_channel(:) = .FALSE.
       flag_ir_channel(:) = .FALSE.

       select case (satellite_id(isatellite))
       case default

          call mpas_log_write("Unrecognized satellite id '"//trim(satellite_id(isatellite))//"'", messageType=MPAS_LOG_CRIT)

       case ("goes18")

          visible_sensor_designation = "v.abi_g18"
          ir_sensor_designation = "abi_g18"

          if (MPAS_field_will_be_written('goes18_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('goes18_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('goes18_brightnessT_channel_08')) flag_ir_channel(8)   = .TRUE.
          if (MPAS_field_will_be_written('goes18_brightnessT_channel_09')) flag_ir_channel(9)   = .TRUE.
          if (MPAS_field_will_be_written('goes18_brightnessT_channel_10')) flag_ir_channel(10)  = .TRUE.

       case ("goes19")

          visible_sensor_designation = "v.abi_g18"  ! I find no v.abi_g19.  Using g18.
          ir_sensor_designation = "abi_g19"

          if (MPAS_field_will_be_written('goes19_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('goes19_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('goes19_brightnessT_channel_08')) flag_ir_channel(8)   = .TRUE.
          if (MPAS_field_will_be_written('goes19_brightnessT_channel_09')) flag_ir_channel(9)   = .TRUE.
          if (MPAS_field_will_be_written('goes19_brightnessT_channel_10')) flag_ir_channel(10)  = .TRUE.

       case ("himawari9")

          visible_sensor_designation = "v.ahi_himawari9"
          ir_sensor_designation = "ahi_himawari9"

          if (MPAS_field_will_be_written('himawari9_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('himawari9_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('himawari9_radiance_channel_03'))    flag_vis_channel(3)  = .TRUE.
          if (MPAS_field_will_be_written('himawari9_brightnessT_channel_08')) flag_ir_channel(8)   = .TRUE.
          if (MPAS_field_will_be_written('himawari9_brightnessT_channel_09')) flag_ir_channel(9)   = .TRUE.
          if (MPAS_field_will_be_written('himawari9_brightnessT_channel_10')) flag_ir_channel(10)  = .TRUE.

       case ("meteosat9")

          visible_sensor_designation = "v.seviri_m09"
          ir_sensor_designation = "seviri_m09"

          if (MPAS_field_will_be_written('meteosat9_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat9_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat9_radiance_channel_12'))    flag_vis_channel(12) = .TRUE.
          if (MPAS_field_will_be_written('meteosat9_brightnessT_channel_05')) flag_ir_channel(5)   = .TRUE.
          if (MPAS_field_will_be_written('meteosat9_brightnessT_channel_06')) flag_ir_channel(6)   = .TRUE.

       case ("meteosat10")

          visible_sensor_designation = "v.seviri_m10"
          ir_sensor_designation = "seviri_m10"

          if (MPAS_field_will_be_written('meteosat10_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat10_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat10_radiance_channel_12'))    flag_vis_channel(12) = .TRUE.
          if (MPAS_field_will_be_written('meteosat10_brightnessT_channel_05')) flag_ir_channel(5)   = .TRUE.
          if (MPAS_field_will_be_written('meteosat10_brightnessT_channel_06')) flag_ir_channel(6)   = .TRUE.

       case ("meteosat11")

          visible_sensor_designation = "v.seviri_m10"  ! I find no visible sensor v.seviri_m11.  Using m10
          ir_sensor_designation = "seviri_m11"

          if (MPAS_field_will_be_written('meteosat11_radiance_channel_01'))    flag_vis_channel(1)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat11_radiance_channel_02'))    flag_vis_channel(2)  = .TRUE.
          if (MPAS_field_will_be_written('meteosat11_radiance_channel_12'))    flag_vis_channel(12) = .TRUE.
          if (MPAS_field_will_be_written('meteosat11_brightnessT_channel_05')) flag_ir_channel(5)   = .TRUE.
          if (MPAS_field_will_be_written('meteosat11_brightnessT_channel_06')) flag_ir_channel(6)   = .TRUE.

       end select

       call setup_sensor_channel_subset(trim(visible_sensor_designation), trim(ir_sensor_designation), &
            flag_vis_channel, flag_ir_channel, max_channels, &
            local_n_sensors, local_sensor_id, &
            sat(isatellite)%channel_subset, sat(isatellite)%nchannels_subset)

       sat(isatellite)%n_sensors = local_n_sensors
       satellite_chinfo_start(isatellite) = n_sensors + 1

       do k = 1, local_n_sensors
          sensor_id(n_sensors + k) = local_sensor_id(k)
       end do
       n_sensors = n_sensors + local_n_sensors

       if (allocated(local_sensor_id)) deallocate(local_sensor_id)

    enddo PRESCAN_SATELLITE_LOOP

    if (n_sensors == 0) then
       return
    endif

    !
    !  Create and fill the atmosphere structure <atm> before any of
    !  our loops for the CRTM work, since we can reuse this.
    !

    allocate(atm(nCells), stat=allocation_status)
    if (allocation_status /= 0) then
       call mpas_log_write("Problem allocating atm structure for CRTM diagnostics", messageType=MPAS_LOG_CRIT)
    endif

    call crtm_atmosphere_create(atm, nVertLevels, n_absorbers, n_clouds, n_aerosols)
    if (any(.not. crtm_atmosphere_associated(atm))) then
       call mpas_log_write("Problem with crtm_atmosphere_create", messageType=MPAS_LOG_CRIT)
    endif

    call fill_crtm_atmosphere_structure(nVertLevels, nCells, atm)

    !
    !  Create and fill the surface structure <sfc> before any of
    !  our loops for the CRTM procedures, since we can reuse this.
    !

    allocate(sfc(nCells), stat=allocation_status)
    if (allocation_status /= 0) then
       call mpas_log_write("Problem allocating sfc structure for CRTM diagnostics", messageType=MPAS_LOG_CRIT)
    endif

    call fill_crtm_surface_structure(nCells, sfc)

    !
    !  Initialize the CRTM coefficients for all the sensors that we
    !  will use subsequently in this routine. chinfo is allocated fresh
    !  here (empty -- nothing to copy in from anywhere) and CRTM_Init
    !  populates it directly, in place. It is never copied again after
    !  this point.
    !

    allocate(chinfo(n_sensors))

    COEFFICIENT_PATH = "./CRTM_NC_COEFFS/"
    IRwaterCoeff_File = "Nalli.IRwater.EmisCoeff.bin"
    IRlandCoeff_File = "IGBP.IRland.EmisCoeff.bin"
    IRsnowCoeff_File = ""
    IRiceCoeff_File = ""
    VISwaterCoeff_File = ""
    VISlandCoeff_File = "IGBP.VISland.EmisCoeff.bin"
    VISsnowCoeff_File = ""
    VISiceCoeff_File = ""
    MWwaterCoeff_File = ""
    ierr = CRTM_Init( sensor_id(1:n_sensors), chinfo, &
         &            CloudCoeff_File = "CloudCoeff.bin", &
         &            File_Path=trim(COEFFICIENT_PATH), &
         &            IRwaterCoeff_File=trim(IRwaterCoeff_File), &
         &            IRlandCoeff_File=trim(IRlandCoeff_File), &
                             !         &            IRsnowCoeff_File=trim(IRsnowCoeff_File), &
                             !         &            IRiceCoeff_File=trim(IRiceCoeff_File), &
                             !         &            VISwaterCoeff_File=trim(VISwaterCoeff_File), &
         &            VISlandCoeff_File=trim(VISlandCoeff_File), &
                             !         &            VISsnowCoeff_File=trim(VISsnowCoeff_File), &
                             !         &            VISiceCoeff_File=trim(VISiceCoeff_File), &
                             !         &            MWwaterCoeff_File=trim(MWwaterCoeff_File), &
         &            Quiet=.TRUE., &
         &            load_cloudcoeff=.TRUE., &
         &            load_aerosolcoeff=.FALSE.)

    if (ierr /= SUCCESS) then
       call mpas_log_write("Problem with CRTM_Init: $i", intArgs=(/ierr/), messageType=MPAS_LOG_CRIT)
    endif

    !
    !  For each satellite and each sensor, subset to the individual
    !  channels we request for each sensor.
    !

    do isatellite = 1, satellite_count
       do n = 1, sat(isatellite)%n_sensors
          idx = satellite_chinfo_start(isatellite) + n - 1
          ierr = crtm_channelinfo_subset(chinfo(idx), &
               channel_subset=sat(isatellite)%channel_subset(n, 1:sat(isatellite)%nchannels_subset(n)))
          if (ierr /= SUCCESS) then
             call mpas_log_write("Problem with CRTM_Channelinfo_subset: $i", intArgs=(/ierr/), messageType=MPAS_LOG_CRIT)
          endif
       enddo
    enddo

    !
    !  Per-satellite work:
    !

    SATELLITE_LOOP : do isatellite = 1, satellite_count

       !
       !  First check: see if we actually need to compute anything for
       !  this satellite.
       !
       
       if (sat(isatellite)%n_sensors == 0) cycle SATELLITE_LOOP

       !
       !  For a given satellite, only a portion of the globe is in
       !  view.  As we fill the <crtm_geometry_type> <geo_subset> in the
       !  routine <fill_crtm_geometry_subset>, we calculate which MPAS
       !  cells are in view of the satellite, and return a
       !  <geo_subset> array that includes only those cells.  The
       !  count of those cells is returned as <nVisible>, and the
       !  indices (of the full global mesh) to that subset of cells is
       !  returned as <visible_idx>.
       !

       !
       !  Create and fill the geometry structure <geo_subset>, returning also <nVisible> and <visible_idx>
       !  to index the cells that are in view of the particular satellite we're working with.
       !

       call fill_crtm_geometry_structure_subset(trim(satellite_id(isatellite)), nCells, &
            &                                   nVisible, visible_idx, geo_subset)

       !
       !  Gather indices that are visible (via <visible_idx>), so we
       !  don't have to make our calculations where the globe is not
       !  in view of the satellite.  <atm_subset> and <sfc_subset> are
       !  now dimensioned by <nVisible>, matching <geo_subset>.
       !

       atm_subset = atm(visible_idx)
       sfc_subset = sfc(visible_idx)

       !
       !  Each satellite we represent may have multiple sensors which
       !  we simulate.  The rest of the work for CRTM is done on a
       !  per-sensor basis in the <SENSOR_LOOP>.  Each pass through
       !  <SENSOR_LOOP> computes CRTM results for all channels
       !  requested for that particular sensor.
       !

       SENSOR_LOOP: do n = 1, sat(isatellite)%n_sensors
          idx = satellite_chinfo_start(isatellite) + n - 1

          n_channels = crtm_channelinfo_n_channels(chinfo(idx))

          !
          !  The <crtm_rtsolution_type> <rts_subset> will hold the
          !  results of the CRTM forward model.  The second dimension
          !  of <rts_subset> we dimension here by <nVisible>, matching
          !  <geo_subset>, <atm_subset>, and <sfc_subset>.  After the
          !  CRTM forward model is called, we distribute the results
          !  from the <nVisible> points in <rts_subset> to the indices
          !  of the full MPAS mesh according to the mapping
          !  represented by <visible_idx>.
          !
          
          !
          !  Set up <rts_subset> to hold the results of the CRTM forward model.
          !
          
          allocate(rts_subset(n_channels, nVisible), stat=allocation_status)
          if (allocation_status /= 0) then
             call mpas_log_write("Problem allocating rts_subset", messageType=MPAS_LOG_CRIT)
          endif
          call mpas_log_write("n_channels for sensor $i = $i", intArgs=(/idx, n_channels/) )

          call crtm_rtsolution_create(rts_subset, nVertLevels)
          if (any(.not. crtm_rtsolution_associated(rts_subset))) then
             call mpas_log_write("Problem with crtm_rtsolution_create", messageType=MPAS_LOG_CRIT)
          endif

          !
          !  Call the CRTM forward model
          !

          call crtm_options_create(opt, n_channels)
          if (ierr /= SUCCESS) then
             call mpas_log_write("Problem with CRTM_options_create: $i", intArgs=(/ierr/), messageType=MPAS_LOG_CRIT)
          endif
          call crtm_options_setvalue(opt, n_Streams=2)
          

          call mpas_log_write("call crtm_forward for sensor $i ("//trim(chinfo(idx)%Sensor_ID)//")", intArgs=(/idx/))

          ierr = crtm_forward(atm_subset, sfc_subset, geo_subset, chinfo(idx:idx), rts_subset, Options=opt)
          if (ierr /= SUCCESS) then
             call mpas_log_write("Problem with CRTM_Forward: $i", intArgs=(/ierr/), messageType=MPAS_LOG_CRIT)
          endif
          call mpas_log_write("Done with CRTM_Forward.")

          !
          !  Distribute CRTM results to the full MPAS diagnostic
          !  arrays.  We build the expected <diagnostic_field_name>
          !  from the satellite name, the sensor channel, and an
          !  expected "radiance" or "brightnessT" category based on
          !  the <sensor_id>.  This latter bit ("radiance" or
          !  "brightnessT") could be a problem: it needs some more
          !  specification since we might down-the-road want to output
          !  radiance (instead of / in addition to) brightnessT from
          !  the IR channels.  TODO.
          !

          if (sensor_id(idx)(1:2) == "v.") then
             diagnostic_field_category = "radiance"
          else
             diagnostic_field_category = "brightnessT"
          endif
          do ichannel = 1, sat(isatellite)%nchannels_subset(n)
             write(diagnostic_field_name, '(A,"_",A,"_channel_",I2.2)') &
                  &     trim(satellite_id(isatellite)), &
                  &     trim(diagnostic_field_category), &
                  &     sat(isatellite)%channel_subset(n,ichannel)
             call mpas_log_write("     Creating field "//trim(diagnostic_field_name))
             call mpas_pool_get_array(diag, trim(diagnostic_field_name), diagnostic_array)
             diagnostic_array = -1.E36
             if (diagnostic_field_category == "radiance") then
                diagnostic_array(visible_idx) = rts_subset(ichannel,:)%radiance
             else if (diagnostic_field_category == "brightnessT") then
                diagnostic_array(visible_idx) = rts_subset(ichannel,:)%brightness_temperature
                             endif
             nullify(diagnostic_array)
             call mpas_log_write("     Done with field "//trim(diagnostic_field_name))
          enddo

          !
          !  Destroy <rts_subset> since we need to recreate it anew for
          !  the next sensor.
          !

          call crtm_rtsolution_destroy(rts_subset)
          deallocate(rts_subset)

       enddo SENSOR_LOOP

       !
       !  Once we are done with all the sensors for a given satellite,
       !  we must destroy the subset arrays that were prepared for that
       !  satellite.
       !

       call crtm_geometry_destroy(geo_subset)
       deallocate(geo_subset)

       call crtm_surface_destroy(sfc_subset)
       deallocate(sfc_subset)

       call crtm_atmosphere_destroy(atm_subset)
       deallocate(atm_subset)

    enddo SATELLITE_LOOP

    !
    !  All done with our computations.  Clean up after ourself,
    !  destroying the <chinfo> structure and the full-mesh <atm> and <sfc>
    !  structures we constructed for this diagnostic output snapshot.
    !

    ierr = crtm_destroy(chinfo)
    if (ierr /= SUCCESS) then
       call mpas_log_write("Problem with CRTM_destroy: $i", intArgs=(/ierr/), messageType=MPAS_LOG_CRIT)
    endif

    call crtm_channelinfo_destroy(chinfo)
    deallocate(chinfo)
    deallocate(sensor_id)

    call crtm_atmosphere_destroy(atm)
    deallocate(atm)

    call crtm_surface_destroy(sfc)
    deallocate(sfc)

    call mpas_log_write("This completes processing in CRTM_DIAGNOSTICS_COMPUTE.")

  end subroutine crtm_diagnostics_compute
  
  !-----------------------------------------------------------------------
  !  routine crtm_diagnostics_reset
  !
  !> \brief Reset diagnostic after it has been written
  !> \author 
  !> \date   
  !> \details
  !>  Reset diagnostic after it has been written
  !
  !-----------------------------------------------------------------------
  subroutine crtm_diagnostics_reset()

    implicit none

    ! Called immediately after diagnostics have been written
    ! Reset counters, accumulations, etc.
    call mpas_log_write("This a call to CRTM_DIAGNOSTICS_RESET, which currently does nothing.")

  end subroutine crtm_diagnostics_reset


  !-----------------------------------------------------------------------
  !  routine crtm_diagnostics_cleanup
  !
  !> \brief Finalizes diagnostic
  !> \author Michael Duda
  !> \date   6 September 2016
  !> \details
  !>  Finalizes diagnostic
  !
  !-----------------------------------------------------------------------
  subroutine crtm_diagnostics_cleanup()

    implicit none

    ! Deallocate scratch arrays, etc.
    call mpas_log_write("This a call to CRTM_DIAGNOSTICS_CLEANUP, which currently does nothing.")

  end subroutine crtm_diagnostics_cleanup

  subroutine setup_sensor_channel_subset(vis_sensor_id, ir_sensor_id, &
       flag_vis_channel, flag_ir_channel, max_channels, &
       n_sensors, sensor_id, channel_subset, nchannels_subset)

    use crtm_module, only: crtm_channelinfo_type

    implicit none

    integer, intent(in) :: max_channels
    
    character(len=*), intent(in) :: vis_sensor_id
    character(len=*), intent(in) :: ir_sensor_id
    logical, dimension(max_channels), intent(in) :: flag_vis_channel
    logical, dimension(max_channels), intent(in) :: flag_ir_channel

    integer, intent(out) :: n_sensors
    character(len=20), allocatable, dimension(:), intent(out) :: sensor_id
    integer, dimension(10,max_channels), intent(inout) :: channel_subset
    integer, dimension(10), intent(inout) :: nchannels_subset

    logical :: use_vis, use_ir
    integer :: isensor, ichannel, allocation_status

    use_vis = any(flag_vis_channel)
    use_ir  = any(flag_ir_channel)

    nchannels_subset(:) = 0
    n_sensors = count([use_vis, use_ir])

    if (n_sensors == 0) return

    allocate(sensor_id(n_sensors), stat=allocation_status)
    if (allocation_status /= 0) then
       call mpas_log_write("Problem allocating structures for CRTM diagnostics", messageType=MPAS_LOG_CRIT)
    endif

    isensor = 0
    if (use_vis) then
       isensor = isensor + 1
       sensor_id(isensor) = vis_sensor_id
       do ichannel = 1, max_channels
          if (flag_vis_channel(ichannel)) then
             nchannels_subset(isensor) = nchannels_subset(isensor) + 1
             channel_subset(isensor, nchannels_subset(isensor)) = ichannel
          endif
       enddo
    endif

    if (use_ir) then
       isensor = isensor + 1
       sensor_id(isensor) = ir_sensor_id
       do ichannel = 1, max_channels
          if (flag_ir_channel(ichannel)) then
             nchannels_subset(isensor) = nchannels_subset(isensor) + 1
             channel_subset(isensor, nchannels_subset(isensor)) = ichannel
          endif
       enddo
    endif

  end subroutine setup_sensor_channel_subset
  
  !-----------------------------------------------------------------------
  !  routine fill_crtm_geometry_structure
  !
  !> \brief fill CRTM geometry structure
  !> \author Kevin W. Manning
  !> \date 10 August 2026
  !> \details
  !>  Utility function for CRTM diagnostics to fill the geometry structure
  !>  with appropriate fields
  !
  !-----------------------------------------------------------------------
  subroutine fill_crtm_geometry_structure(satellite_id, nCells, geo)
    use mpas_pool_routines, only : mpas_pool_get_array
    use mpas_pool_routines, only : mpas_pool_get_config
    use mpas_timekeeping  , only : MPAS_NOW
    use mpas_timekeeping  , only : mpas_get_clock_time
    use mpas_timekeeping  , only : mpas_get_time
    use crtm_module       , only : crtm_geometry_type
    implicit none
    character(len=*), intent(in) :: satellite_id
    integer, intent(in) :: nCells
    type(crtm_geometry_type)   , dimension(nCells) :: geo

    integer :: ierr
    real (kind=RKIND) :: r_sat
    real (kind=RKIND) :: sensor_zenith_angle
    real (kind=RKIND) :: sensor_azimuth_angle
    real (kind=RKIND) :: cos_gamma
    real (kind=RKIND) :: d
    real (kind=RKIND) :: x
    real (kind=RKIND) :: y
    real (kind=RKIND) :: dlon
    real (kind=RKIND) :: gamma
    real (kind=RKIND) :: rlat0
    real (kind=RKIND) :: rlon0
    real (kind=RKIND) :: satellite_elevation_meters
    real (kind=RKIND) :: satellite_latitude
    real (kind=RKIND) :: satellite_longitude
    real (kind=RKIND) :: scan_angle
    real (kind=RKIND) :: sin_scan_angle
    real (kind=RKIND) :: cos_zenith_angle
    type (MPAS_time_type) :: currTime
    integer :: mpas_day_of_year
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    integer :: second

    character(len=StrKIND), pointer :: xtime
    real (kind=RKIND)     , pointer :: sphere_radius
    real (kind=RKIND)     , dimension(:)    , pointer :: latCell
    real (kind=RKIND)     , dimension(:)    , pointer :: lonCell
    real (kind=RKIND)     , dimension(:,:)  , pointer :: zgrid
    real (kind=RKIND)     , dimension(:)    , pointer :: coszr
    real (kind=RKIND) :: cos_sz
    real (kind=RKIND) :: doy
    real (kind=RKIND) :: gamma_time
    real (kind=RKIND) :: sol_decl
    real (kind=RKIND) :: eq_time
    real (kind=RKIND) :: utc_fractional_hours
    real (kind=RKIND) :: az_sun
    
    integer :: iCell

    call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
    call mpas_pool_get_array(mesh, 'latCell', latCell)
    call mpas_pool_get_array(mesh, 'lonCell', lonCell)
    call mpas_pool_get_array(mesh, 'zgrid', zgrid)
    call mpas_pool_get_array(diag_physics, 'coszr' ,coszr )

    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write("Problem with call to mpas_get_clock_time from CRTM diagnostics module", messageType=MPAS_LOG_CRIT)
    endif
    call mpas_get_time(currTime, YYYY=year, MM=month, DD=day, H=hour, M=minute, S=second, DoY=mpas_day_of_year, ierr=ierr)
    utc_fractional_hours = real(hour, RKIND) + real(minute, RKIND)/60.0 + real(second, RKIND)/3600.0
    doy = real(mpas_day_of_year, kind=RKIND) + utc_fractional_hours/24.0

    ! call mpas_pool_get_array(state, 'xtime', xtime)
    ! read(xtime(1:4), '(i4)') year
    ! read(xtime(6:7), '(i2)') month
    ! read(xtime(9:10), '(i2)') day
    ! read(xtime(12:13), '(i2)') hour
    ! read(xtime(15:16), '(i2)') minute
    ! read(xtime(18:19), '(i2)') second
    call mpas_log_write("CRTM_DIAGNOSTICS for $i, $i, $i, $i, $i, $i, $r", &
         &              intArgs=(/year, month, day, hour, minute, second/), &
         &              realArgs=(/doy/))
    
    
    !
    !  Fill the geometry structure 'geo'
    !

    if (satellite_id == "goes18") then
       satellite_latitude = 0.0
       satellite_longitude = 223.0 ! 137.0 West Longitude
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    elseif (satellite_id == "goes19") then
       satellite_latitude = 0.0
       satellite_longitude = 284.8 ! 75.2 West Longitude
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    elseif (satellite_id == "himawari9") then
       satellite_latitude = 0.0
       satellite_longitude = 140.7
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    elseif (satellite_id == "meteosat9") then
       satellite_latitude = 0.0
       satellite_longitude = 45.5
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    elseif (satellite_id == "meteosat10") then
       satellite_latitude = 0.0
       satellite_longitude = 0.0
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    elseif (satellite_id == "meteosat11") then
       satellite_latitude = 0.0
       satellite_longitude = 9.5
       satellite_elevation_meters = 1.E3 * 35786  ! meters above the sphere
    else
       call mpas_log_write("Unrecognized satellite: '"//satellite_id//"'", messageType=MPAS_LOG_CRIT)
    endif

    rlat0 = satellite_latitude * degrad
    rlon0 = satellite_longitude * degrad
    r_sat = sphere_radius + satellite_elevation_meters

    if (.TRUE.) then

       gamma_time = (2.0 * pii / 365.0) * (doy - 1.0 + (utc_fractional_hours - 12.0) / 24.0)

       !
       !  Solar Declination Angle (in radians)
       !
       sol_decl = 0.006918 - 0.399912*cos(gamma_time) + 0.070257*sin(gamma_time) &
            &   - 0.006758*cos(2.0*gamma_time) + 0.000907*sin(2.0*gamma_time) &
            &   - 0.002697*cos(3.0*gamma_time) + 0.001480*sin(3.0*gamma_time)

       !
       !  Equation of Time (in minutes)
       !
       eq_time = 229.18 * ( 0.000075 + 0.001868*cos(gamma_time) - 0.032077*sin(gamma_time) &
            &    - 0.014615*cos(2.0*gamma_time) - 0.040849*sin(2.0*gamma_time) )
    endif



    do iCell = 1, nCells
       geo(iCell)%latitude = latCell(iCell) / degrad
       geo(iCell)%longitude = lonCell(iCell) / degrad
       if (geo(iCell)%longitude < 0) then
          geo(iCell)%longitude = geo(iCell)%longitude + 360
       endif

       cos_gamma = sin(rlat0)*sin(latCell(iCell)) + cos(rlat0)*cos(latCell(iCell))*cos(lonCell(iCell)-rlon0)
       cos_gamma = MAX(-1.0, MIN(1.0, cos_gamma))
       gamma = acos(cos_gamma)

       d = sqrt( &
            &    sphere_radius**2 &
            &    + r_sat**2 &
            &    - 2 * sphere_radius * r_sat * cos_gamma )
       sin_scan_angle = sphere_radius*sin(gamma)/d
       sin_scan_angle = MAX(-1.0, MIN(1.0, sin_scan_angle))
       scan_angle = asin(sin_scan_angle) / degrad


       cos_zenith_angle = (d**2 + sphere_radius**2 - r_sat**2) / (2.0 * d * sphere_radius)
       cos_zenith_angle = MAX(-1.0, MIN(1.0, cos_zenith_angle))
       sensor_zenith_angle = (pii - acos(cos_zenith_angle)) / degrad

       if ( (abs(rlat0 - latCell(iCell)) < 1.E-8) .and. (abs(rlon0 - lonCell(iCell)) < 1.E-8)) then
          sensor_azimuth_angle = 0.0
       else
          dlon = rlon0-lonCell(iCell)
          y = sin(dlon)*cos(rlat0)
          x = (cos(latCell(iCell))*sin(rlat0)-sin(latCell(iCell))*cos(rlat0)*cos(dlon))
          sensor_azimuth_angle = atan2(y,x)/degrad
          if (sensor_azimuth_angle<0.0) sensor_azimuth_angle = sensor_azimuth_angle + 360.0
       endif

       if (sensor_zenith_angle > 89.99) then
          sensor_zenith_angle = 0.0
          sensor_azimuth_angle = 999.9
       else if (sensor_zenith_angle > 80) then
          sensor_zenith_angle = 80.0
          sin_scan_angle = sphere_radius * sin(sensor_zenith_angle*degrad)/r_sat
          sin_scan_angle = MAX(-1.0, MIN(1.0, sin_scan_angle))
          scan_angle = asin(sin_scan_angle) / degrad
       endif

       if (.TRUE.) then
          block
            real(kind=RKIND) :: lon_deg
            real(kind=RKIND) :: sin_sz
            real(kind=RKIND) :: cos_sa
            real(kind=RKIND) :: hr_angle
            real(kind=RKIND) :: local_solar_time
            
            lon_deg = lonCell(iCell) / degrad
            if (lon_deg > 180.0) lon_deg = lon_deg - 360.0 ! Normalize to -180 to +180

            !  Local Solar Time in minutes
            local_solar_time = utc_fractional_hours * 60.0 + eq_time + 4.0 * lon_deg

            !  Hour angle in radians (0 at local solar noon, negative in morning, positive in afternoon)
            hr_angle = ((local_solar_time / 4.0) - 180.0) * degrad

            !
            !  Spherical Trig for Solar Azimuth
            !
            cos_sz = sin(latCell(iCell)) * sin(sol_decl) + cos(latCell(iCell)) * cos(sol_decl) * cos(hr_angle)
            cos_sz = MAX(-1.0, MIN(1.0, cos_sz))
            sin_sz = sqrt(1.0 - cos_sz**2)

            if (sin_sz < 1.E-6) then
               az_sun = 0.0 ! Sun is directly at the zenith (rare, polar, or sub-solar point exact match)
            else
               ! Calculate solar azimuth using atan2 to properly resolve quadrants
               y = -sin(hr_angle) * cos(sol_decl)
               x = cos(latCell(iCell)) * sin(sol_decl) - sin(latCell(iCell)) * cos(sol_decl) * cos(hr_angle)
               az_sun = atan2(y, x) / degrad

               ! Map from [-180, 180] to CRTM convention [0, 360] degrees clockwise from North
               if (az_sun < 0.0) az_sun = az_sun + 360.0
            endif
          end block
       else
          az_sun = 0.0
       endif

       geo(iCell)%surface_altitude = zgrid(1,iCell)

       ! geo(iCell)%sensor_scan_angle = 0.0     ! Assume each MPAS cell has a satellite sensor directly overhead
       ! geo(iCell)%sensor_zenith_angle = 0.0   ! Assume each MPAS cell has a satellite sensor directly overhead
       ! geo(iCell)%sensor_azimuth_angle = 0.0  ! Assume each MPAS cell has a satellite sensor directly overhead
       ! geo(iCell)%source_zenith_angle = 0.0   ! Assume each MPAS cell has the sun directly overhead

       geo(iCell)%sensor_scan_angle = scan_angle
       geo(iCell)%sensor_zenith_angle = sensor_zenith_angle
       geo(iCell)%sensor_azimuth_angle = sensor_azimuth_angle
       ! geo(iCell)%source_zenith_angle = ACOS( MAX(-1.0, MIN(1.0, coszr(iCell))) ) * 180.0 / pii
       geo(iCell)%source_zenith_angle = ACOS( MAX(-1.0, MIN(1.0, cos_sz)) ) * 180.0 / pii
       geo(iCell)%source_azimuth_angle = az_sun

       geo(iCell)%iFOV = 0 ! Default value.  Used for things like antenna correction.

       geo(iCell)%year = year
       geo(iCell)%month = month
       geo(iCell)%day = day
    enddo

  end subroutine fill_crtm_geometry_structure

!-----------------------------------------------------------------------
  !  routine fill_crtm_geometry_structure_subset
  !
  !> \brief fill CRTM geometry structure, visible-cell subset only
  !> \details
  !>  Same formulas as the original fill_crtm_geometry_structure, but
  !>  restructured so the expensive back half of the per-cell computation
  !>  (sensor azimuth angle, solar geometry block) only runs for cells
  !>  that pass the cheap visibility test, and geo_subset/visible_idx are
  !>  built directly -- no full nCells-sized geo array is ever allocated
  !>  or filled.
  !-----------------------------------------------------------------------
  subroutine fill_crtm_geometry_structure_subset(satellite_id, nCells, &
       nVisible, visible_idx, geo_subset)
    use mpas_pool_routines, only : mpas_pool_get_array
    use mpas_pool_routines, only : mpas_pool_get_config
    use mpas_timekeeping  , only : MPAS_NOW
    use mpas_timekeeping  , only : mpas_get_clock_time
    use mpas_timekeeping  , only : mpas_get_time
    use crtm_module       , only : crtm_geometry_type
    implicit none
    character(len=*), intent(in) :: satellite_id
    integer, intent(in) :: nCells
    integer, intent(out) :: nVisible
    integer, allocatable, dimension(:), intent(out) :: visible_idx
    type(crtm_geometry_type), allocatable, dimension(:), intent(out) :: geo_subset

    integer :: ierr
    real (kind=RKIND) :: r_sat
    real (kind=RKIND) :: sensor_zenith_angle
    real (kind=RKIND) :: sensor_azimuth_angle
    real (kind=RKIND) :: cos_gamma
    real (kind=RKIND) :: d
    real (kind=RKIND) :: x
    real (kind=RKIND) :: y
    real (kind=RKIND) :: dlon
    real (kind=RKIND) :: gamma
    real (kind=RKIND) :: rlat0
    real (kind=RKIND) :: rlon0
    real (kind=RKIND) :: satellite_elevation_meters
    real (kind=RKIND) :: satellite_latitude
    real (kind=RKIND) :: satellite_longitude
    real (kind=RKIND) :: scan_angle
    real (kind=RKIND) :: sin_scan_angle
    real (kind=RKIND) :: cos_zenith_angle
    type (MPAS_time_type) :: currTime
    integer :: mpas_day_of_year
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    integer :: second

    character(len=StrKIND), pointer :: xtime
    real (kind=RKIND)     , pointer :: sphere_radius
    real (kind=RKIND)     , dimension(:)    , pointer :: latCell
    real (kind=RKIND)     , dimension(:)    , pointer :: lonCell
    real (kind=RKIND)     , dimension(:,:)  , pointer :: zgrid
    real (kind=RKIND)     , dimension(:)    , pointer :: coszr
    real (kind=RKIND) :: cos_gamma_threshold
    real (kind=RKIND) :: cos_rlat0
    real (kind=RKIND) :: sin_rlat0
    real (kind=RKIND) :: cos_sz
    real (kind=RKIND) :: doy
    real (kind=RKIND) :: gamma_time
    real (kind=RKIND) :: sol_decl
    real (kind=RKIND) :: eq_time
    real (kind=RKIND) :: utc_fractional_hours
    real (kind=RKIND) :: az_sun

    logical, allocatable, dimension(:) :: bitmask
    integer :: iCell, i, allocation_status

    call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
    call mpas_pool_get_array(mesh, 'latCell', latCell)
    call mpas_pool_get_array(mesh, 'lonCell', lonCell)
    call mpas_pool_get_array(mesh, 'zgrid', zgrid)
    call mpas_pool_get_array(diag_physics, 'coszr' ,coszr )

    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write("Problem with call to mpas_get_clock_time from CRTM diagnostics module", messageType=MPAS_LOG_CRIT)
    endif
    call mpas_get_time(currTime, YYYY=year, MM=month, DD=day, H=hour, M=minute, S=second, DoY=mpas_day_of_year, ierr=ierr)
    utc_fractional_hours = real(hour, RKIND) + real(minute, RKIND)/60.0 + real(second, RKIND)/3600.0
    doy = real(mpas_day_of_year, kind=RKIND) + utc_fractional_hours/24.0

    call mpas_log_write("CRTM_DIAGNOSTICS for $i, $i, $i, $i, $i, $i, $r", &
         &              intArgs=(/year, month, day, hour, minute, second/), &
         &              realArgs=(/doy/))

    !
    !  Satellite position, unchanged from the original
    !

    if (satellite_id == "goes18") then
       satellite_latitude = 0.0
       satellite_longitude = 223.0 ! 137.0 West Longitude
       satellite_elevation_meters = 1.E3 * 35786
    elseif (satellite_id == "goes19") then
       satellite_latitude = 0.0
       satellite_longitude = 284.8 ! 75.2 West Longitude
       satellite_elevation_meters = 1.E3 * 35786
    elseif (satellite_id == "himawari9") then
       satellite_latitude = 0.0
       satellite_longitude = 140.7
       satellite_elevation_meters = 1.E3 * 35786
    elseif (satellite_id == "meteosat9") then
       satellite_latitude = 0.0
       satellite_longitude = 45.5
       satellite_elevation_meters = 1.E3 * 35786
    elseif (satellite_id == "meteosat10") then
       satellite_latitude = 0.0
       satellite_longitude = 0.0
       satellite_elevation_meters = 1.E3 * 35786
    elseif (satellite_id == "meteosat11") then
       satellite_latitude = 0.0
       satellite_longitude = 9.5
       satellite_elevation_meters = 1.E3 * 35786
    else
       call mpas_log_write("Unrecognized satellite: '"//satellite_id//"'", messageType=MPAS_LOG_CRIT)
    endif

    rlat0 = satellite_latitude * degrad
    rlon0 = satellite_longitude * degrad
    sin_rlat0 = sin(rlat0)
    cos_rlat0 = cos(rlat0)
    r_sat = sphere_radius + satellite_elevation_meters

    gamma_time = (2.0 * pii / 365.0) * (doy - 1.0 + (utc_fractional_hours - 12.0) / 24.0)

    sol_decl = 0.006918 - 0.399912*cos(gamma_time) + 0.070257*sin(gamma_time) &
         &   - 0.006758*cos(2.0*gamma_time) + 0.000907*sin(2.0*gamma_time) &
         &   - 0.002697*cos(3.0*gamma_time) + 0.001480*sin(3.0*gamma_time)

    eq_time = 229.18 * ( 0.000075 + 0.001868*cos(gamma_time) - 0.032077*sin(gamma_time) &
         &    - 0.014615*cos(2.0*gamma_time) - 0.040849*sin(2.0*gamma_time) )

    !
    !  Precompute once, before the loop: the cos_gamma visibility
    !  threshold corresponding to the 89.99-degree zenith cutoff, derived
    !  by inverting the zenith-angle formula below at Z=89.99 deg
    !
    
    block
      real(kind=RKIND) :: Z_thresh, cosZ, sinZ, d_thresh
      Z_thresh = 89.99 * degrad
      cosZ = cos(Z_thresh)
      sinZ = sin(Z_thresh)
      d_thresh = -sphere_radius*cosZ + sqrt(r_sat**2 - sphere_radius**2 * sinZ**2)
      cos_gamma_threshold = (sphere_radius**2 + r_sat**2 - d_thresh**2) / (2.0*sphere_radius*r_sat)
    end block

    !
    !  Pass 1: cheap visibility test only -- cos_gamma/gamma/d/sensor_zenith_angle,
    !  no azimuth, no solar geometry. Runs over all nCells, but this is the
    !  minimum work needed to know which cells are visible at all.
    !

    allocate(bitmask(nCells))
    do iCell = 1, nCells
       cos_gamma = sin_rlat0*sin(latCell(iCell)) + cos_rlat0*cos(latCell(iCell))*cos(lonCell(iCell)-rlon0)
       bitmask(iCell) = (cos_gamma >= cos_gamma_threshold)
    enddo

    visible_idx = pack( [(i, i=1,nCells)], bitmask )
    nVisible = size(visible_idx)
    deallocate(bitmask)

    !
    !  Pass 2: full geometry, only for the nVisible cells that passed the
    !  visibility test. This recomputes the cheap chain above (small cost,
    !  now over nVisible cells instead of nCells) and then does the
    !  expensive azimuth + solar geometry work that pass 1 skipped
    !  entirely for invisible cells.
    !

    allocate(geo_subset(nVisible), stat=allocation_status)
    if (allocation_status /= 0) then
       call mpas_log_write("Problem allocating geo_subset for CRTM diagnostics", messageType=MPAS_LOG_CRIT)
    endif

    do i = 1, nVisible
       iCell = visible_idx(i)

       geo_subset(i)%latitude = latCell(iCell) / degrad
       geo_subset(i)%longitude = lonCell(iCell) / degrad
       if (geo_subset(i)%longitude < 0) then
          geo_subset(i)%longitude = geo_subset(i)%longitude + 360
       endif

       cos_gamma = sin_rlat0*sin(latCell(iCell)) + cos_rlat0*cos(latCell(iCell))*cos(lonCell(iCell)-rlon0)
       cos_gamma = MAX(-1.0, MIN(1.0, cos_gamma))
       gamma = acos(cos_gamma)

       d = sqrt( sphere_radius**2 + r_sat**2 - 2 * sphere_radius * r_sat * cos_gamma )
       sin_scan_angle = sphere_radius*sin(gamma)/d
       sin_scan_angle = MAX(-1.0, MIN(1.0, sin_scan_angle))
       scan_angle = asin(sin_scan_angle) / degrad

       cos_zenith_angle = (d**2 + sphere_radius**2 - r_sat**2) / (2.0 * d * sphere_radius)
       cos_zenith_angle = MAX(-1.0, MIN(1.0, cos_zenith_angle))
       sensor_zenith_angle = (pii - acos(cos_zenith_angle)) / degrad

       ! sensor_zenith_angle <= 89.99 is guaranteed here (that's why iCell
       ! is in visible_idx), so the ">89.99 -> invisible sentinel" branch
       ! from the original never applies in pass 2 -- only the >80 clamp
       ! can still fire.

       if ( (abs(rlat0 - latCell(iCell)) < 1.E-8) .and. (abs(rlon0 - lonCell(iCell)) < 1.E-8)) then
          sensor_azimuth_angle = 0.0
       else
          dlon = rlon0-lonCell(iCell)
          y = sin(dlon)*cos_rlat0
          x = (cos(latCell(iCell))*sin_rlat0-sin(latCell(iCell))*cos_rlat0*cos(dlon))
          sensor_azimuth_angle = atan2(y,x)/degrad
          if (sensor_azimuth_angle<0.0) sensor_azimuth_angle = sensor_azimuth_angle + 360.0
       endif

       if (sensor_zenith_angle > 80) then
          sensor_zenith_angle = 80.0
          sin_scan_angle = sphere_radius * sin(sensor_zenith_angle*degrad)/r_sat
          sin_scan_angle = MAX(-1.0, MIN(1.0, sin_scan_angle))
          scan_angle = asin(sin_scan_angle) / degrad
       endif

       block
         real(kind=RKIND) :: lon_deg
         real(kind=RKIND) :: sin_sz
         real(kind=RKIND) :: cos_sa
         real(kind=RKIND) :: hr_angle
         real(kind=RKIND) :: local_solar_time

         lon_deg = lonCell(iCell) / degrad
         if (lon_deg > 180.0) lon_deg = lon_deg - 360.0

         local_solar_time = utc_fractional_hours * 60.0 + eq_time + 4.0 * lon_deg
         hr_angle = ((local_solar_time / 4.0) - 180.0) * degrad

         cos_sz = sin(latCell(iCell)) * sin(sol_decl) + cos(latCell(iCell)) * cos(sol_decl) * cos(hr_angle)
         cos_sz = MAX(-1.0, MIN(1.0, cos_sz))
         sin_sz = sqrt(1.0 - cos_sz**2)

         if (sin_sz < 1.E-6) then
            az_sun = 0.0
         else
            y = -sin(hr_angle) * cos(sol_decl)
            x = cos(latCell(iCell)) * sin(sol_decl) - sin(latCell(iCell)) * cos(sol_decl) * cos(hr_angle)
            az_sun = atan2(y, x) / degrad
            if (az_sun < 0.0) az_sun = az_sun + 360.0
         endif
       end block

       geo_subset(i)%surface_altitude = zgrid(1,iCell)
       geo_subset(i)%sensor_scan_angle = scan_angle
       geo_subset(i)%sensor_zenith_angle = sensor_zenith_angle
       geo_subset(i)%sensor_azimuth_angle = sensor_azimuth_angle
       geo_subset(i)%source_zenith_angle = ACOS( MAX(-1.0, MIN(1.0, cos_sz)) ) * 180.0 / pii
       geo_subset(i)%source_azimuth_angle = az_sun
       geo_subset(i)%iFOV = 0

       geo_subset(i)%year = year
       geo_subset(i)%month = month
       geo_subset(i)%day = day
    enddo

  end subroutine fill_crtm_geometry_structure_subset  
  
  !-----------------------------------------------------------------------
  !  routine fill_crtm_surface_structure
  !
  !> \brief fill CRTM surface structure
  !> \author Kevin W. Manning
  !> \date 10 August 2026
  !> \details
  !>  Utility function for CRTM diagnostics to fill the surface structure
  !>  with appropriate fields
  !
  !-----------------------------------------------------------------------
  subroutine fill_crtm_surface_structure(nCells, sfc)
    use mpas_pool_routines, only : mpas_pool_get_array
    use crtm_module, only: crtm_surface_type
    implicit None
    integer, intent(in) :: nCells
    type(crtm_surface_type), intent(inout)   , dimension(nCells) :: sfc
    
    integer, pointer :: iswater
    integer, pointer :: isice
    integer,           dimension(:)    , pointer :: ivgtyp
    real(kind=RKIND),  dimension(:)    , pointer :: xicem
    ! real(kind=RKIND),  dimension(:)    , pointer :: snowc
    integer,           dimension(:)    , pointer :: landmask
    real(kind=RKIND),  dimension(:)    , pointer :: smstot
    ! real(kind=RKIND),  dimension(:)    , pointer :: canwat
    real(kind=RKIND),  dimension(:)    , pointer :: vegfra
    real(kind=RKIND),  dimension(:)    , pointer :: skintemp
    real(kind=RKIND),  dimension(:)    , pointer :: lai
    real(kind=RKIND),  dimension(:)    , pointer :: snowh
    real(kind=RKIND),  dimension(:)    , pointer :: u10
    real(kind=RKIND),  dimension(:)    , pointer :: v10
    real(kind=RKIND),  dimension(:,:)  , pointer :: tslb
    integer :: iCell

    ! My guess for mapping MPAS soil categories to CRTM soil type classification.
    ! integer, dimension(19), parameter :: remap_soil_type = (/ 1, 1, 4, 2, 2, 6, 7, 6, 6, 5, 3, 3, 8, 1, 1, 9, 1, 1, 1 /)

    
    call mpas_pool_get_array(sfc_input, 'iswater', iswater)
    call mpas_pool_get_array(sfc_input, 'isice', isice)
    call mpas_pool_get_array(sfc_input, 'ivgtyp', ivgtyp)
    call mpas_pool_get_array(diag_physics, 'xicem', xicem )
    call mpas_pool_get_array(sfc_input, 'landmask', landmask )
    call mpas_pool_get_array(diag_physics, 'smstot', smstot )
    ! call mpas_pool_get_array(diag_physics, 'canwat', canwat )
    call mpas_pool_get_array(sfc_input, 'skintemp', skintemp )
    call mpas_pool_get_array(sfc_input, 'vegfra', vegfra )
    call mpas_pool_get_array(diag_physics, 'lai', lai )
    call mpas_pool_get_array(sfc_input, 'tslb', tslb )
    call mpas_pool_get_array(sfc_input, 'snowh', snowh )
    call mpas_pool_get_array(diag_physics, 'u10', u10 )
    call mpas_pool_get_array(diag_physics, 'v10', v10 )

    if ((iswater == 17) .and. (isice == 15)) then
       ! modified igbp
    else
       call mpas_log_write("For crtm, iswater = $i and isice = $i", intArgs=(/iswater, isice/), messageType=MPAS_LOG_CRIT)
    endif

    !
    !  Fill the surface structure
    !

    do iCell = 1, nCells
       if ((ivgtyp(iCell) > 20).or.(ivgtyp(iCell)<1)) then
          call mpas_log_write("For crtm, ivgtype = $i is a problem", intArgs=(/ivgtyp(iCell)/), messageType=MPAS_LOG_CRIT)
       endif
       if (landmask(iCell) == 1) then
          ! Land point (or permanent ice shelf)
          if (ivgtyp(iCell) == iswater) then
             call mpas_log_write("For crtm, iswater at land", messageType=MPAS_LOG_CRIT)
          endif
          
          if (ivgtyp(iCell) == isice) then
             ! Call it land with "old snow"?
             sfc(iCell)%land_coverage = 0.0
             sfc(iCell)%snow_coverage = 1.0
             sfc(iCell)%water_coverage = 0.0
             sfc(iCell)%ice_coverage = 0.0
             sfc(iCell)%snow_type = 1 ! 1 == "old snow"
          else
             ! land point.  We may want to look at the snowc field.
             ! if (snowc(iCell) > 0) then
             !    sfc(iCell)%snow_coverage = 1.0
             !    sfc(iCell)%land_coverage = 0.0
             ! else
             !    sfc(iCell)%snow_coverage = 0.0
             !    sfc(iCell)%land_coverage = 1.0
             ! endif
             sfc(iCell)%snow_coverage = 0.0
             sfc(iCell)%land_coverage = 1.0
             sfc(iCell)%water_coverage = 0.0
             sfc(iCell)%ice_coverage = 0.0
             sfc(iCell)%snow_type = 2 ! 2 == "new snow"
          endif
          
       else if (landmask(iCell) == 0) then
          ! Water/Sea-ice point
          sfc(iCell)%land_coverage = 0.0
          sfc(iCell)%water_coverage = 1.0
          sfc(iCell)%ice_coverage = 0.0
          
          if (xicem(iCell) > 0.0) then
             sfc(iCell)%water_coverage = 1.0-xicem(iCell)
             sfc(iCell)%ice_coverage = xicem(iCell)
          else
             sfc(iCell)%water_coverage = 1.0
             sfc(iCell)%ice_coverage = 0.0
          endif
             
          sfc(iCell)%snow_coverage = 0.0
          sfc(iCell)%water_type = 1
          sfc(iCell)%snow_type = 2
          sfc(iCell)%ice_type = 1
          sfc(iCell)%land_type = iswater

       endif

       if (abs(1.0-(sfc(iCell)%land_coverage + sfc(iCell)%ice_coverage + sfc(iCell)%snow_coverage + sfc(iCell)%water_coverage)) > 1.E-6) then
          call mpas_log_write("landmask(iCell) = $i", intArgs=(/landmask(iCell)/))
          call mpas_log_write("For crtm, $i coverage: $r $r $r $r", intArgs = (/landmask(iCell)/), &
               & realArgs=(/real(sfc(iCell)%land_coverage, kind=4), real(sfc(iCell)%ice_coverage, kind=4), real(sfc(iCell)%snow_coverage, kind=4), real(sfc(iCell)%water_coverage, kind=4)/), messageType=MPAS_LOG_CRIT)
       endif


       sfc(iCell)%land_temperature = skintemp(iCell)
       sfc(iCell)%water_temperature = skintemp(iCell)
       sfc(iCell)%snow_temperature = skintemp(iCell)
       sfc(iCell)%ice_temperature = skintemp(iCell)
       sfc(iCell)%soil_temperature = tslb(1,iCell)

       
       sfc(iCell)%soil_moisture_content = smstot(iCell) ! smstot in m3/m3 is the same as g/cm3

       ! Need more information to convert MPAS canwat (kg/m2) to CRTM canopy_water_content (g/cm3)
       ! So for now, just use the CRTM default value
       sfc(iCell)%canopy_water_content = 0.05 ! canwat(iCell)
       sfc(iCell)%vegetation_fraction = vegfra(iCell)
       sfc(iCell)%lai = lai(iCell)

       ! For microwave channels only?
       ! sfc(iCell)%soil_type = remap_soil_type(isltyp(iCell))

       sfc(iCell)%soil_type = 8
       
       sfc(iCell)%wind_speed = hypot(u10(iCell), v10(iCell))
       sfc(iCell)%wind_direction = atan2(-u10(iCell), -v10(iCell)) / degrad
       if (sfc(iCell)%wind_direction < 0) then
          sfc(iCell)%wind_direction = sfc(iCell)%wind_direction + 360
       endif
       
       ! sfc(iCell)%salinity =
       sfc(iCell)%snow_depth = snowh(iCell) * 1.E3  ! convert m to mm
       ! sfc(iCell)%snow_density = 
       ! sfc(iCell)%snow_grain_size = 
       ! sfc(iCell)%ice_thickness =
       ! sfc(iCell)%ice_density = 
       ! sfc(iCell)%ice_roughness = 
       ! sfc(iCell)%SensorData = 
          
    enddo
    
  end subroutine fill_crtm_surface_structure

  !-----------------------------------------------------------------------
  !  routine fill_crtm_atmosphere_structure
  !
  !> \brief fill CRTM atmosphere structure
  !> \author Kevin W. Manning
  !> \date 10 August 2026
  !> \details
  !>  Utility function for CRTM diagnostics to fill the atmosphere structure
  !>  with appropriate fields
  !
  !-----------------------------------------------------------------------
  subroutine fill_crtm_atmosphere_structure(nVertLevels, nCells, atm)
    use mpas_pool_routines, only : mpas_pool_get_dimension
    use mpas_pool_routines, only : mpas_pool_get_array
    use mpas_timekeeping  , only : MPAS_NOW
    use mpas_timekeeping  , only : mpas_get_clock_time
    use mpas_timekeeping  , only : mpas_get_time
    use crtm_module, only: crtm_atmosphere_type

    ! CRTM "Climatology" options
    use crtm_module, only: TROPICAL
    use crtm_module, only: MIDLATITUDE_SUMMER
    use crtm_module, only: MIDLATITUDE_WINTER
    use crtm_module, only: SUBARCTIC_SUMMER
    use crtm_module, only: SUBARCTIC_WINTER
    use crtm_module, only: US_STANDARD_ATMOSPHERE
    use crtm_module, only: MASS_MIXING_RATIO_UNITS
    use crtm_module, only: VOLUME_MIXING_RATIO_UNITS
    use crtm_module, only: WATER_CLOUD
    use crtm_module, only: ICE_CLOUD
    use crtm_module, only: SNOW_CLOUD
    use crtm_module, only: RAIN_CLOUD
    use crtm_module, only: GRAUPEL_CLOUD
    use crtm_module, only: H2O_ID
    use crtm_module, only: O3_ID
    implicit None
    integer, intent(in) :: nVertLevels
    integer, intent(in) :: nCells
    type(crtm_atmosphere_type), intent(inout)   , dimension(nCells) :: atm
    type (MPAS_time_type) :: currTime

    integer, pointer :: index_qv
    integer, pointer :: index_qc
    integer, pointer :: index_qr
    integer, pointer :: index_qi
    integer, pointer :: index_qs
    integer, pointer :: index_qg
    integer, pointer :: index_ni
    integer, pointer :: index_nr
    integer, pointer :: index_nc
    real (kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=RKIND), dimension(:,:)  , pointer :: cldfrac
    real (kind=RKIND), dimension(:,:)  , pointer :: zgrid
    real (kind=RKIND), dimension(:,:)  , pointer :: theta_m
    real (kind=RKIND), dimension(:,:)  , pointer :: exner
    real (kind=RKIND), dimension(:,:)  , pointer :: pressure_p
    real (kind=RKIND), dimension(:,:)  , pointer :: pressure_base
    real (kind=RKIND), dimension(:,:)  , pointer :: o3vmr
    ! real (kind=RKIND), dimension(:,:)  , pointer :: re_cloud
    ! real (kind=RKIND), dimension(:,:)  , pointer :: re_ice
    ! real (kind=RKIND), dimension(:,:)  , pointer :: re_snow
    real (kind=RKIND), dimension(:)    , pointer :: latCell
    ! real (kind=RKIND), dimension(:)    , pointer :: lonCell

    real (kind=RKIND), allocatable, dimension(:) :: p8w
    real (kind=RKIND), allocatable, dimension(:) :: p
    real (kind=RKIND) :: qsum
    real (kind=RKIND) :: q
    real (kind=RKIND) :: delta_p
    real (kind=RKIND) :: z0
    real (kind=RKIND) :: z1
    real (kind=RKIND) :: z2
    real (kind=RKIND) :: w1
    real (kind=RKIND) :: w2

    integer :: ierr
    integer :: iCell
    integer :: k
    integer :: season
    integer :: mpas_day_of_year
    real (kind=RKIND) :: dlat
    integer :: imp_physics

    call mpas_pool_get_array(mesh, 'latCell', latCell)
    ! call mpas_pool_get_array(mesh, 'lonCell', lonCell)
    call mpas_pool_get_dimension(state, 'index_qv', index_qv)
    call mpas_pool_get_dimension(state, 'index_qc', index_qc)
    call mpas_pool_get_dimension(state, 'index_qr', index_qr)
    call mpas_pool_get_dimension(state, 'index_qi', index_qi)
    call mpas_pool_get_dimension(state, 'index_qs', index_qs)
    call mpas_pool_get_dimension(state, 'index_qg', index_qg)
    call mpas_pool_get_dimension(state, 'index_ni', index_ni)
    call mpas_pool_get_dimension(state, 'index_nr', index_nr)
    call mpas_pool_get_dimension(state, 'index_nc', index_nc)
    call mpas_pool_get_array(state, 'scalars', scalars, 1)
    call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
    call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
    call mpas_pool_get_array(diag, 'exner', exner)
    call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
    call mpas_pool_get_array(mesh, 'zgrid', zgrid)
    call mpas_pool_get_array(diag_physics, 'o3vmr', o3vmr)
    call mpas_pool_get_array(diag_physics, 'cldfrac', cldfrac)
    ! call mpas_pool_get_array(diag_physics, 're_cloud' ,re_cloud )
    ! call mpas_pool_get_array(diag_physics, 're_ice' ,re_ice )
    ! call mpas_pool_get_array(diag_physics, 're_snow' ,re_snow )

    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write("Problem with call to mpas_get_clock_time from CRTM diagnostics module", messageType=MPAS_LOG_CRIT)
    endif
    call mpas_get_time(currTime, doy=mpas_day_of_year, ierr=ierr)

    if ((mpas_day_of_year > 79) .and. (mpas_day_of_year < 266)) then
       season = 1  ! N Hemisphere Summer, S Hemisphere winter;  with spring and fall equinoctes around day 79 and 266
    else
       season = -1 ! N Hemisphere Winter, S Hemisphere summer
    endif


    allocate(p8w(nVertLevels+1))
    allocate(p(nVertLevels))

    do iCell = 1, nCells

       do k = 1, nVertLevels
          p(k) = (pressure_base(k,iCell) + pressure_p(k,iCell))
       enddo

       !calculation of total pressure at cell centers (at vertical velocity points):
       k = nVertLevels+1
       z0 = zgrid(k,iCell)
       z1 = 0.5*(zgrid(k,iCell)+zgrid(k-1,iCell)) 
       z2 = 0.5*(zgrid(k-1,iCell)+zgrid(k-2,iCell))
       w1 = (z0-z2)/(z1-z2)
       w2 = 1.-w1
       !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
       p8w(k) = exp(w1*log(p(k-1))+w2*log(p(k-2)))

       do k = 2, nVertLevels
          w1 = (zgrid(k,iCell)-zgrid(k-1,iCell)) / (zgrid(k+1,iCell)-zgrid(k-1,iCell))
          w2 = (zgrid(k+1,iCell)-zgrid(k,iCell)) / (zgrid(k+1,iCell)-zgrid(k-1,iCell))
          p8w(k) = exp(w1*log(p(k))+w2*log(p(k-1)))
       enddo

       k = 1
       z0 = zgrid(k,iCell)
       z1 = 0.5*(zgrid(k,iCell)+zgrid(k+1,iCell)) 
       z2 = 0.5*(zgrid(k+1,iCell)+zgrid(k+2,iCell))
       w1 = (z0-z2)/(z1-z2)
       w2 = 1.-w1
       p8w(k) = exp(w1*log(p(k))+w2*log(p(k+1)))

       ! atm(iCell)%climatology = US_STANDARD_ATMOSPHERE

       dlat = season * latCell(iCell) / degrad
       if (abs(dlat) < 30) then
          atm(iCell)%climatology = TROPICAL
       else
          if (dlat > 0) then
             ! Latitude times Season will be positive in summer, for either N or S hemisphere.
             if (dlat < 55) then
                atm(iCell)%climatology = MIDLATITUDE_SUMMER
             else
                atm(iCell)%climatology = SUBARCTIC_SUMMER
             endif
          else
             ! Latitude times Season will be negative in winter, for either N or S hemisphere.
             if (dlat > -55) then
                atm(iCell)%climatology = MIDLATITUDE_WINTER
             else
                atm(iCell)%climatology = SUBARCTIC_WINTER
             endif
          endif
       endif

       atm(iCell)%absorber_id = (/H2O_ID, O3_ID/)
       atm(iCell)%absorber_units = (/MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS/)

       ! CRTM expects top-down.
       do k = 1, nVertLevels
          atm(iCell)%pressure(nVertlevels+1-k) = 1.E-2 * p(k)
       enddo

       do k = 1, nVertLevels+1
          atm(iCell)%level_pressure(nVertLevels+2-k) = 1.E-2 * p8w(k)
       enddo

       ! temperature
       do k = 1, nVertLevels
          atm(iCell)%temperature(nVertLevels+1-k) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell)
       enddo

       ! H2O and ozone
       do k = 1, nVertLevels
          atm(iCell)%absorber(nVertLevels+1-k,1) = 1.E3*scalars(index_qv,k,iCell)! H2O values; 1.E3 factor converts from MPAS kg/kg to CRTM g/kg
          if (atm(iCell)%absorber(nVertLevels+1-k,1) < 1.E-15) then
             atm(iCell)%absorber(nVertLevels+1-k,1) = 1.E-15
          endif
          atm(iCell)%absorber(nVertLevels+1-k,2) = o3vmr(k,iCell) * 1.E6
       enddo

       !
       !  Clouds
       !
       
       if (trim(config_microp_scheme) == "mp_wsm6") then
          imp_physics = 6
          atm(iCell)%cloud(1)%type = WATER_CLOUD
          atm(iCell)%cloud(2)%type = ICE_CLOUD
          atm(iCell)%cloud(3)%type = SNOW_CLOUD
          atm(iCell)%cloud(4)%type = RAIN_CLOUD
          atm(iCell)%cloud(5)%type = GRAUPEL_CLOUD
          do k = 1, nVertLevels
             delta_p = p8w(k)-p8w(k+1)
             atm(iCell)%cloud_fraction(nVertLevels+1-k) = cldfrac(k,iCell)

             qsum = scalars(index_qv,k,iCell) + &
                  & scalars(index_qc,k,iCell) + &
                  & scalars(index_qr,k,iCell) + &
                  & scalars(index_qi,k,iCell) + &
                  & scalars(index_qs,k,iCell) + &
                  & scalars(index_qg,k,iCell)
             q = scalars(index_qv,k,iCell)/(1.+qsum)

             !
             !  WATER_CLOUD
             !

             atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) = scalars(index_qc,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(1)%effective_radius(nVertLevels+1-k) = re_cloud(k,iCell)*1.E6          
             atm(iCell)%cloud(1)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  imp_physics,'C')
             if (atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) = 1.E-25
             endif

             !
             !  ICE_CLOUD
             !

             atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) = scalars(index_qi,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(2)%effective_radius(nVertLevels+1-k) = re_ice(k,iCell)*1.E6
             atm(iCell)%cloud(2)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  imp_physics,'I')

             if (atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) = 1.E-25
             endif

             !
             !  SNOW_CLOUD
             !

             atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) = scalars(index_qs,k,iCell)*delta_p/9.80665
             atm(iCell)%cloud(3)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  imp_physics,'S')
             ! atm(iCell)%cloud(3)%effective_radius(nVertLevels+1-k) = re_snow(k,iCell)*1.E6
             if (atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) = 1.E-25
             endif


             !
             !  RAIN_CLOUD
             !

             atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) = scalars(index_qr,k,iCell)*delta_p/9.80665
             atm(iCell)%cloud(4)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  imp_physics,'R')
             ! atm(iCell)%cloud(4)%effective_radius(nVertLevels+1-k) = re_rain(k,iCell)*1.E6
             if (atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) = 1.E-25
             endif


             !
             !  GRAUPEL_CLOUD
             !

             atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) = scalars(index_qg,k,iCell)*delta_p/9.80665
             atm(iCell)%cloud(5)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  -1.E36, &
                  -1.E36, &
                  -1.E36, &
                  imp_physics,'G')
             ! atm(iCell)%cloud(5)%effective_radius(nVertLevels+1-k) = re_graupel(k,iCell)*1.E6
             if (atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) = 1.E-25
             endif
          enddo
       else if (trim(config_microp_scheme) == "mp_thompson") then
          imp_physics = 8
          atm(iCell)%cloud(1)%type = WATER_CLOUD
          atm(iCell)%cloud(2)%type = ICE_CLOUD
          atm(iCell)%cloud(3)%type = SNOW_CLOUD
          atm(iCell)%cloud(4)%type = RAIN_CLOUD
          atm(iCell)%cloud(5)%type = GRAUPEL_CLOUD
          do k = 1, nVertLevels
             delta_p = p8w(k)-p8w(k+1)
             atm(iCell)%cloud_fraction(nVertLevels+1-k) = cldfrac(k,iCell)

             qsum = scalars(index_qv,k,iCell) + &
                  & scalars(index_qc,k,iCell) + &
                  & scalars(index_qr,k,iCell) + &
                  & scalars(index_qi,k,iCell) + &
                  & scalars(index_qs,k,iCell) + &
                  & scalars(index_qg,k,iCell)
             q = scalars(index_qv,k,iCell)/(1.+qsum)

             !
             !  WATER CLOUD
             !
             
             atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) = scalars(index_qc,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(1)%effective_radius(nVertLevels+1-k) = re_cloud(k,iCell)*1.E6          
             atm(iCell)%cloud(1)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &  ! f_rimef unused for Thompson
                  -1.E36, &  ! nlice unused for Thompson
                  -1.E36, &  ! nrain unused for Thompson
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  scalars(index_nr,k,iCell), &
                  scalars(index_ni,k,iCell), &
                  scalars(index_nc,k,iCell), &
                  imp_physics,'C')
             if (atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(1)%water_content(nVertLevels+1-k) = 1.E-25
             endif
             
             !
             !  ICE CLOUD
             !
             
             atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) = scalars(index_qi,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(2)%effective_radius(nVertLevels+1-k) = re_ice(k,iCell)*1.E6          
             atm(iCell)%cloud(2)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &  ! f_rimef unused for Thompson
                  -1.E36, &  ! nlice unused for Thompson
                  -1.E36, &  ! nrain unused for Thompson
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  scalars(index_nr,k,iCell), &
                  scalars(index_ni,k,iCell), &
                  scalars(index_nc,k,iCell), &
                  imp_physics,'I')
             if (atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(2)%water_content(nVertLevels+1-k) = 1.E-25
             endif

             !
             !  SNOW CLOUD
             !
             
             atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) = scalars(index_qs,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(3)%effective_radius(nVertLevels+1-k) = re_snow(k,iCell)*1.E6          
             atm(iCell)%cloud(3)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &  ! f_rimef unused for Thompson
                  -1.E36, &  ! nlice unused for Thompson
                  -1.E36, &  ! nrain unused for Thompson
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  scalars(index_nr,k,iCell), &
                  scalars(index_ni,k,iCell), &
                  scalars(index_nc,k,iCell), &
                  imp_physics,'S')
             if (atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(3)%water_content(nVertLevels+1-k) = 1.E-25
             endif

             !
             !  RAIN CLOUD
             !
             
             atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) = scalars(index_qr,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(4)%effective_radius(nVertLevels+1-k) = re_rain(k,iCell)*1.E6          
             atm(iCell)%cloud(4)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &  ! f_rimef unused for Thompson
                  -1.E36, &  ! nlice unused for Thompson
                  -1.E36, &  ! nrain unused for Thompson
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  scalars(index_nr,k,iCell), &
                  scalars(index_ni,k,iCell), &
                  scalars(index_nc,k,iCell), &
                  imp_physics,'R')
             if (atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(4)%water_content(nVertLevels+1-k) = 1.E-25
             endif

             !
             !  GRAUPEL CLOUD
             !
             
             atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) = scalars(index_qg,k,iCell)*delta_p/9.80665
             ! atm(iCell)%cloud(5)%effective_radius(nVertLevels+1-k) = re_graupel(k,iCell)*1.E6          
             atm(iCell)%cloud(5)%effective_radius(nVertLevels+1-k) = effr(&
                  p(k), &
                  real(atm(iCell)%temperature(nVertLevels+1-k), kind=4), &
                  q, &
                  scalars(index_qc,k,iCell), &
                  scalars(index_qi,k,iCell), &
                  scalars(index_qr,k,iCell), &
                  -1.E36, &  ! f_rimef unused for Thompson
                  -1.E36, &  ! nlice unused for Thompson
                  -1.E36, &  ! nrain unused for Thompson
                  scalars(index_qs,k,iCell), &
                  scalars(index_qg,k,iCell), &
                  scalars(index_nr,k,iCell), &
                  scalars(index_ni,k,iCell), &
                  scalars(index_nc,k,iCell), &
                  imp_physics,'G')
             if (atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) < 1.E-25) then
                atm(iCell)%cloud(5)%water_content(nVertLevels+1-k) = 1.E-25
             endif
          enddo
       endif

       ! I have no aerosols for now.

    enddo
  end subroutine fill_crtm_atmosphere_structure

  !-------------------------------------------------------------------------------
!> @brief EFFR Computes effective particle radii channel selection using LVLS from WRF_CNTRL.PARM. 
!>
!> @param pmid real Mid-layer pressure.
!> @param t real Temperature.
!> @param q real Specific humidity.
!> @param qqw real Cloud water mixing ratio.
!> @param qqi real Ice mixing ratio.
!> @param qqr real Rain mixing ratio.
!> @param f_rimef real "Rime Factor", ratio of total ice growth to deposition growth.
!> @param nlice real Time-averaged number concentration of large ice.
!> @param nrain real Number concentration of rain drops.
!> @param qqs real Snow mixing ratio.
!> @param qqg real Graupel mixing ratio.
!> @param qqnr real Rain number concentration.
!> @param qqni real Ice number concentration.
!> @param qqnw real cloud water number concentration.
!> @param mp_opt integer Microphysics option.
!> @param species character Particle type (e.g., cloud, rain, graupel, snow, ice).
!> @return EFFR Effective particle radii channel selection. 
!-------------------------------------------------------------------------------

  REAL FUNCTION EFFR(pmid,t,q,qqw,qqi,qqr,f_rimef, nlice, nrain, &
       qqs,qqg,qqnr,qqni,qqnw,mp_opt,species)

    !       JASON OTKIN AND WILLIAM LEWIS
    !       09 DECEMBER 2014
    !       Greg Thompson, 20200924
    !
    !       Imported from UPP

    ! use params_mod, only: pi, rd, d608, rg

    implicit none

    real, parameter :: PI=3.141592653589793     !< Pi
    real, parameter :: RD=287.04                !< Gas constant for dry air
    real, parameter :: D608=0.608               !< Decimal number 0.608
    real, parameter :: G=9.81                   !< Acceleration due to gravity
    real, parameter :: GI=1./G                  !< Inverse of acceleration due to gravity
    real, parameter :: RG=1./GI                 !< Reciprocal of GI (equal to G)
    

    real :: pmid,t,q,qqw,qqi,qqr,qqs,qqg,f_rimef,nlice,nrain
    real :: qqnr,qqni,qqnw
    character(LEN=1) :: species

    integer                         :: n,count,count1,mp_opt
    real :: rho, ncc, rhox
    real :: n0_s, n0_r, n0_g
    real :: lambdar, lambdas, lambdag

    !-------------------------------------------------------------------------------
    !  GAMMA FUNCTION & RELATED VARIABLES
    !-------------------------------------------------------------------------------

    real :: gamma
    real :: gamma_crg, gamma_s
    !       real :: gamma_i

    ! real :: WGAMMA, GAMMLN

    real    :: rc,am_c,bm_c,cce(3,15),ccg(3,15),ocg1(15),ocg2(15)
    integer :: nu_c

    real, dimension(0:15), parameter:: g_ratio = (/6,24,60,120,210, &
         &              336,504,720,990,1320,1716,2184,2730,3360,4080,4896/)

    real    :: rr, mu_r, am_r, bm_r, cre(3), crg(3), ore1, org1, org2
    real    :: mvd_r, ron_sl, ron_r0, ron_c0, ron_c1, ron_c2, obmr

    real    :: ri, mu_i, am_i, bm_i, cie(3), cig(3), oig1, oig2, obmi

    real    :: rs, am_s, oams, cse(3)
    real    :: loga, a, b, tc0, smob, smo2, smoc
    REAL, PARAMETER:: mu_s = 0.6357
    REAL, PARAMETER:: Kap0 = 490.6
    REAL, PARAMETER:: Kap1 = 17.46
    REAL, PARAMETER:: Lam0 = 20.78
    REAL, PARAMETER:: Lam1 = 3.29

    !-------------------------------------------------------------------------------
    !  MINIMUM/MAXIMUM CONSTANTS FOR ALL SCHEMES
    !-------------------------------------------------------------------------------

    real, parameter :: eps=0.622, beta_crg=3., beta_i=2.,beta_s=2.4

    real, parameter :: min_qc=1.e-7, min_qr=1.e-7, min_qi=1.e-8,min_qs=1.e-8, min_qg=1.e-7
    real, parameter :: min_c=2.e-6,  min_r=20.e-6, min_i=4.e-6,min_s=20.e-6, min_g=20.e-6
    real, parameter :: max_c=1.e-2,  max_r=1.e-2,  max_i=1.e-3,max_s=2.e-2,  max_g=5.e-0

    real    :: am_g, bm_g, mu_g
    real    :: cgg(3), cge(3), oge1, obmg, ogg1, ogg2

    real    :: ygra1, zans1, rg2
    double precision :: no_exp, no_min, lm_exp, lamg, lamc, lamr, lami, lams

    !-------------------------------------------------------------------------------
    !  WSM6-SPECIFIC ARRAYS
    !-------------------------------------------------------------------------------

    real :: wsm6_nci, xmi, xmitemp

    !-------------------------------------------------------------------------------
    !  CONSTANTS FOR WSM6 MICROPHYSICS SCHEME
    !-------------------------------------------------------------------------------

    real, parameter :: wsm6_cnp=3.e8, wsm6_rhor=1000.
    real, parameter :: wsm6_rhos=100., wsm6_rhog=500.
    real, parameter :: wsm6_dimax=500.e-6, wsm6_dicon=11.9
    real, parameter :: wsm6_alpha=.12, wsm6_t0c=273.15
    real, parameter :: wsm6_n0s=2.e6, wsm6_n0smax=1.e11
    real, parameter :: wsm6_n0r=8.e6, wsm6_n0g=4.e6
    real, parameter :: wsm6_qmin=1.e-15

    !-------------------------------------------------------------------------------
    !  CONSTANTS FOR LIN MICROPHYSICS SCHEME
    !-------------------------------------------------------------------------------

    real, parameter :: lin_rhoi=100., lin_rhor=1000., lin_rhos=100.
    real, parameter :: lin_rhog=400., lin_cnp=3.e8
    real, parameter :: lin_n0r=8.e6,  lin_n0s=3.e6,   lin_n0g=4.e6

    !-------------------------------------------------------------------------------
    !  CONSTANTS FOR NEW THOMPSON MICROPHYSICS SCHEME (FOR WRF VERSIONS 3.1 AND UP)
    !-------------------------------------------------------------------------------

    real, parameter :: nthom_rhor=1000., nthom_rhos=100.
    !       WM LEWIS updated rhog to 500 from 400
    real, parameter :: nthom_rhog=500.,  nthom_rhoi=890.
    real, parameter :: nthom_gon_min=1.e4, nthom_gon_max=3.e6
    real, parameter :: nthom_nt_c=100.e6

    real, parameter :: nthom_min_nci=5.e2
    real, parameter :: nthom_min_ncr=1.e-6

    real, parameter :: nthom_bm_s=2.0               !this is important

    real :: nci2, ncc2, ncr2

    real, dimension(10), parameter :: &
         nthom_sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
         0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
    real, dimension(10), parameter :: &
         nthom_sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
         0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

    !-------------------------------------------------------------------------------
    !  CONSTANTS FOR GFDL MICROPHYSICS SCHEME - which is Lin for precip clouds
    !-------------------------------------------------------------------------------

    real, parameter :: gfdl_rhoi=100., gfdl_rhor=1000., gfdl_rhos=100.
    real, parameter :: gfdl_rhog=400., gfdl_cnp=3.e8
    real, parameter :: gfdl_tice = 273.16

    real, parameter :: gfdl_qmin = 1.0e-5, gfdl_ccn = 1.0e8, gfdl_beta = 1.22
    real, parameter :: gfdl_gammar = 17.837789, gfdl_gammas = 8.2850630, gfdl_gammag = 11.631769
    real, parameter :: gfdl_alphar = 0.8, gfdl_alphas = 0.25, gfdl_alphag = 0.5
    real, parameter :: gfdl_n0r=8.e6, gfdl_n0s=3.e6, gfdl_n0g=4.e6

    real, parameter :: gfdl_rewmin = 5.0,  gfdl_rewmax = 10.0
    real, parameter :: gfdl_reimin = 10.0, gfdl_reimax = 150.0
    real, parameter :: gfdl_rermin = 0.0,  gfdl_rermax = 10000.0
    real, parameter :: gfdl_resmin = 0.0,  gfdl_resmax = 10000.0
    real, parameter :: gfdl_regmin = 0.0,  gfdl_regmax = 10000.0



    if(mp_opt==6) then                        !WSM6 SCHEME

       n0_r = wsm6_n0r
       n0_g = wsm6_n0g
       n0_s = wsm6_n0s

    elseif(mp_opt==2)then                     !LIN SCHEME

       n0_r = lin_n0r
       n0_g = lin_n0g
       n0_s = lin_n0s

    endif

    gamma_crg = 6.0      ! gamma(1.0 + beta_crg)
    gamma_s = 2.981134   ! gamma(1.0 + beta_s)
    !       gamma_i = 2.0        ! gamma(1.0 + beta_i)

    !------------------------------------------------------------------------------
    !  SET DIAMETER ARRAYS TO ZERO, COMPUTE DENSITY
    !------------------------------------------------------------------------------

    effr=0. 

    rho=pmid/(rd*t*(1.+D608*q))


    if(mp_opt==6)then

       SELECT CASE(species)

       CASE("C")

          if ( qqw>min_qc ) then !cloud diameter: assume constant # concentration
             effr = 1.0E6*(( 6. * rho * qqw ) / &
                  (pi * wsm6_rhor * wsm6_cnp))**(1/3.)

          endif

       CASE("R")

          if ( qqr>min_qr ) then !rain diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqr ) / &
                  ( pi * wsm6_rhor * n0_r * gamma_crg ) ) ** (1/(1+beta_crg ) )
          endif

       CASE("G")

          if ( qqg>min_qg ) then !graupel diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqg ) / &
                  ( pi * wsm6_rhog * n0_g * gamma_crg ) ) ** (1/(1+beta_crg ) )
          endif

       CASE("S")

          if ( qqs>min_qs ) then !snow diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqs ) / &
                  ( pi * wsm6_rhos * n0_s * gamma_s   ) ) ** ( 1/(1+beta_s) )
          endif

          !  ICE DIAMETER: CALCULATED USING METHOD OUTLINED IN WRF BROWSER.  Refer to
          !  phys/module_mp_wsm6.F (Vice:fallout of ice crystal).

       CASE("I")

          if ( qqi>min_qi ) then !ice diameter
             !       wsm6_nci = min(max(5.38e7*(rho*max(qqi,wsm6_qmin)),1.e3),1.e6)
             !       xmi = rho * qqi / wsm6_nci
             !       effr = 1.0E6*min( sqrt(xmi), wsm6_dimax)
             !!      from wsm6, HWRF ver 3.6:
             !!          temp = (den(i,k)*max(qci(i,k,2),qmin))
             !!          temp = sqrt(sqrt(temp*temp*temp))
             !!          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
             !!      diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
             xmitemp=rho*max(qqi,wsm6_qmin)
             xmitemp=sqrt(sqrt(xmitemp*xmitemp*xmitemp))
             xmi= min(max(5.38e7*xmitemp,1.e3),1.e6)
             effr = 1.0E6*max(min(wsm6_dicon * sqrt(xmi),wsm6_dimax), 1.e-25)
          endif

       END SELECT

    elseif(mp_opt==2)then

       SELECT CASE(species)

       CASE("C")

          if ( qqw > min_qc ) then !cloud diameter: assume constant # concentration
             effr = 1.0E6*(( 6. * rho * qqw ) / &
                  (pi * lin_rhor * lin_cnp))**(1/3.)
          endif

       CASE("R")

          if ( qqr > min_qr ) then !rain diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqr ) / &
                  ( pi * lin_rhor * n0_r * gamma_crg ) ) ** (1/(1+beta_crg ) )
          endif

       CASE("I")

          if ( qqi > min_qi ) then !ice diameter: assume constant # concentrtion
             effr = 1.0E6*( ( 6. * rho * qqi ) / &
                  ( pi * lin_rhoi * lin_cnp ) ) ** ( 1/3.)
          endif

       CASE("S")

          if ( qqs > min_qs ) then !snow diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqs ) / &
                  ( pi * lin_rhos * n0_s * gamma_s   ) ) ** ( 1/(1+beta_s) )
          endif

       CASE("G")

          if ( qqg > min_qg ) then !graupel diameter: assume gamma distribution
             effr = 1.0E6*( ( 6. * rho * qqg ) / &
                  ( pi * lin_rhog * n0_g * gamma_crg ) ) ** (1/(1+beta_crg ) )
          endif

       END SELECT

    elseif(mp_opt==8 .or. mp_opt==28)then

       !  rain section

       bm_r   = 3.0
       mu_r   = 0.0
       obmr   = 1.0 / bm_r
       am_r   = pi * nthom_rhor / 6.0

       cre(1) = bm_r + 1.
       cre(2) = mu_r + 1.
       cre(3) = bm_r + mu_r + 1.

       crg(1) = WGAMMA(cre(1))
       crg(2) = WGAMMA(cre(2))
       crg(3) = WGAMMA(cre(3))

       ore1   = 1. / cre(1)
       org1   = 1. / crg(1)
       org2   = 1. / crg(2)

       !  cloud section

       bm_c   = bm_r

       do n = 1, 15
          cce(1,n) = n + 1.             ! Substitute variable value of mu_c
          cce(2,n) = bm_c + n + 1.      ! Substitute variable value of mu_c

          ccg(1,n) = WGAMMA(cce(1,n))
          ccg(2,n) = WGAMMA(cce(2,n))

          ocg1(n)   = 1./ccg(1,n)
          ocg2(n)   = 1./ccg(2,n)
       enddo

       !  ice section

       am_i   = pi * nthom_rhoi / 6.0
       bm_i   = 3.0
       mu_i   = 0.

       cie(1) = mu_i + 1.
       cie(2) = bm_i + mu_i + 1.

       cig(1) = WGAMMA(cie(1))
       cig(2) = WGAMMA(cie(2))

       oig1   = 1./cig(1)
       oig2   = 1./cig(2)
       obmi   = 1./bm_i

       !  snow section

       am_s   = 0.069

       oams   = 1./am_s

       cse(1) = nthom_bm_s + 1.

       !  graupel section
       bm_g   = 3.0
       mu_g   = 0.0
       obmg   = 1.0 / bm_g
       am_g   = pi * nthom_rhog / 6.0

       cge(1) = bm_g + 1.
       cge(2) = mu_g + 1.
       cge(3) = bm_g + mu_g + 1.

       cgg(1) = WGAMMA(cge(1))
       cgg(2) = WGAMMA(cge(2))
       cgg(3) = WGAMMA(cge(3))

       oge1   = 1. / cge(1)
       ogg1   = 1. / cgg(1)
       ogg2   = 1. / cgg(2)

       !CLOUD DIAMETER CALCULATIONS

       SELECT CASE (species)

       CASE("C")

          if(qqw >= min_qc) then

             rc = MAX(1.E-12, qqw * rho)

             if (mp_opt==8) then
                ncc2 = nthom_nt_c
             elseif (mp_opt==28) then
                ncc2 = MAX(1.E-6, qqnw * rho)
             endif

             if (ncc2 < 10.e6) then
                nu_c = 15
             else
                nu_c   = min (15, NINT(1000.e6/ncc2) + 2)
             endif

             lamc = (ncc2/rc)**obmr * (am_r*g_ratio(nu_c))**obmr

             effr = 1.0E6*MAX(4.01E-6, MIN(SNGL(1.0D0*DBLE(3.+nu_c)/lamc),50.E-6))

             !           old UPP
             !             effr = 2.*10.

          endif

          !RAIN DIAMETER CALCULATIONS

       CASE("R")

          if( qqr > min_qr) then

             rr = MAX(1.E-12, qqr * rho)
             ncr2 = MAX(1.E-6, qqnr * rho)
             lamr = (ncr2/rr)**obmr * (am_r*crg(3)*org2)**obmr

             effr = 1.0E6*MAX(50.01E-6, MIN(SNGL(1.0D0*DBLE(3.+mu_r)/lamr),1999.E-6))

             !             old UPP
             !              effr=2.*200.

             !              print*,'effr_rain=',effr/2.

          endif

          !ICE DIAMETER CACLULATIONS

       CASE("I")

          if(qqi >= min_qi) then

             ri = MAX(1.E-12, qqi * rho)
             nci2 = MAX(1.E-6, qqni * rho)

             lami = (nci2/ri)**obmi * (am_i*cig(2)*oig1)**obmi

             effr = 1.0E6*MAX(10.01E-6, MIN(SNGL(1.0D0*DBLE(3.+mu_i)/lami),250.E-6))

             !             old UPP
             !               effr=2.*25.

          endif

          !SNOW DIAMETER CALCULATIONS

       CASE("S")

          rs = MAX(1.E-12, qqs * rho)

          if(qqs >= min_qs) then

             tc0  = min(-0.1, t-273.15)
             smob = rs*oams

             if (nthom_bm_s>(2.0-1.e-3) .and. nthom_bm_s<(2.0+1.e-3))then
                smo2 = smob
             else
                loga = nthom_sa(1) + nthom_sa(2)*tc0 + nthom_sa(3)*nthom_bm_s+               &
                     nthom_sa(4)*tc0*nthom_bm_s + nthom_sa(5)*tc0*tc0 +&
                     nthom_sa(6)*nthom_bm_s*nthom_bm_s +nthom_sa(7)*tc0*tc0*nthom_bm_s +   &
                     nthom_sa(8)*tc0*nthom_bm_s*nthom_bm_s +nthom_sa(9)*tc0*tc0*tc0 +      &
                     nthom_sa(10)*nthom_bm_s*nthom_bm_s*nthom_bm_s

                a    = 10.0**loga

                b    = nthom_sb(1) + nthom_sb(2)*tc0 + nthom_sb(3)*nthom_bm_s+               &
                     nthom_sb(4)*tc0*nthom_bm_s + nthom_sb(5)*tc0*tc0 +&
                     nthom_sb(6)*nthom_bm_s*nthom_bm_s +nthom_sb(7)*tc0*tc0*nthom_bm_s +   &
                     nthom_sb(8)*tc0*nthom_bm_s*nthom_bm_s +nthom_sb(9)*tc0*tc0*tc0 +      &
                     nthom_sb(10)*nthom_bm_s*nthom_bm_s*nthom_bm_s
                smo2 = (smob/a)**(1./b)
             endif

             !Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
             loga      = nthom_sa(1) + nthom_sa(2)*tc0 + nthom_sa(3)*cse(1) +&
                  nthom_sa(4)*tc0*cse(1) + nthom_sa(5)*tc0*tc0 +&
                  nthom_sa(6)*cse(1)*cse(1) + nthom_sa(7)*tc0*tc0*cse(1)+      &
                  nthom_sa(8)*tc0*cse(1)*cse(1) +nthom_sa(9)*tc0*tc0*tc0 +     &
                  nthom_sa(10)*cse(1)*cse(1)*cse(1)

             a       = 10.0**loga

             b       = nthom_sb(1)+ nthom_sb(2)*tc0 + nthom_sb(3)*cse(1) +&
                  nthom_sb(4)*tc0*cse(1) + nthom_sb(5)*tc0*tc0 +&
                  nthom_sb(6)*cse(1)*cse(1) + nthom_sb(7)*tc0*tc0*cse(1) +&
                  nthom_sb(8)*tc0*cse(1)*cse(1) + nthom_sb(9)*tc0*tc0*tc0 +       &
                  nthom_sb(10)*cse(1)*cse(1)*cse(1)

             smoc      = a * smo2**b

             effr = 1.0E6*MAX(50.E-6, MIN(smoc/smob, 1999.E-6))

             !              print*,'snow effr=',effr

             !             changing snow effr recovers "old" UPP Thompson almost exactly;
             !             i.e. the snow effr is the source of the cold discprepancy.

             !             old UPP
             !              effr=2.*250.

          endif

       CASE("G")

          if(qqg >= min_qg) then

             rg2 = MAX(1.E-12, qqg * rho)

             ygra1 = alog10(max(1.E-9, rg2))

             zans1 = 3. + 2./7. * (ygra1+7.)
             zans1 = MAX(2., MIN(zans1, 7.))

             no_exp = 10.**(zans1)

             lm_exp = (no_exp*am_g*cgg(1)/rg2)**oge1

             lamg = lm_exp * (cgg(3)*ogg2*ogg1)**obmg

             effr= 1.0E6*MAX(99.E-6, MIN(SNGL((3.0+mu_g)/lamg), 9999.E-6))

             !           old UPP
             !            effr=350.

          endif

       END SELECT

    elseif(mp_opt==11)then ! GFDL 

       SELECT CASE(species)

       CASE("C")

          ! cloud water (martin et al., 1994)
          if (qqw > min_qc) then
             effr = exp (1.0 / 3.0 * log ((3. * qqw ) / (4. * pi * gfdl_rhor * gfdl_ccn))) * 1.0e6
             effr = max (gfdl_rewmin, min (gfdl_rewmax, effr))
             effr = effr*2.  ! because need diameter here, converted to radius at exit
          end if

       CASE("I")

          ! cloud ice (heymsfield and mcfarquhar, 1996)
          if (qqi > min_qi) then
             if ((t-gfdl_tice) < - 50) then
                effr = gfdl_beta / 9.917 * exp ((1 - 0.891) * log (1.0e3 * qqi)) * 1.0e3
             elseif ((t-gfdl_tice) < - 40.) then
                effr = gfdl_beta / 9.337 * exp ((1 - 0.920) * log (1.0e3 * qqi)) * 1.0e3
             elseif ((t-gfdl_tice) < - 30.) then
                effr = gfdl_beta / 9.208 * exp ((1 - 0.945) * log (1.0e3 * qqi)) * 1.0e3
             else
                effr = gfdl_beta / 9.387 * exp ((1 - 0.969) * log (1.0e3 * qqi)) * 1.0e3
             endif
             effr = max (gfdl_reimin, min (gfdl_reimax, effr))
             effr = effr*2.  ! because need diameter here, converted to radius at exit
          end if

       CASE("R")

          if ( qqr > min_qr ) then !rain diameter: assume gamma distribution
             lambdar = exp (0.25 * log (pi * gfdl_rhor * gfdl_n0r / qqr))
             effr  =  0.5*exp (log (gfdl_gammar / 6.) / gfdl_alphar) / lambdar * 1.0e6
             effr = max (gfdl_rermin, min (gfdl_rermax, effr))
             effr = effr*2.  ! because need diameter here, converted to radius at exit
          endif


       CASE("S")

          if ( qqs > min_qs ) then !snow diameter: assume gamma distribution
             lambdas = exp (0.25 * log (pi * gfdl_rhos * gfdl_n0s / qqs))
             effr = 0.5 * exp (log (gfdl_gammas / 6.) / gfdl_alphas) / lambdas * 1.0e6
             effr = max (gfdl_resmin, min (gfdl_resmax, effr))
             effr = effr*2.  ! because need diameter here, converted to radius at exit
          endif

       CASE("G")

          if ( qqg > min_qg ) then !graupel diameter: assume gamma distribution
             lambdag = exp (0.25 * log (pi * gfdl_rhog * gfdl_n0g / qqg))
             effr = 0.5 * exp (log (gfdl_gammag / 6.) / gfdl_alphag) / lambdag * 1.0e6
             effr = max (gfdl_regmin, min (gfdl_regmax, effr))
             effr = effr*2.  ! because need diameter here, converted to radius at exit
          endif

       END SELECT


    elseif(mp_opt==5.or.mp_opt==85.or.mp_opt==95)then

       SELECT CASE (species)

       CASE("C")

          effr=2.*10.

       CASE("I")

          effr=2.*25.

       CASE("R")

          if( qqr > min_qr) then
             rhox=1000.
             effr=2.*1.0E6*1.5*(rho*qqr/(pi*rhox*nrain))**(1./3.)

             !      old UPP    
             !      effr=2.*200.
             !      effr=min(200.,effr)
             !      print*,'effr_rain=',effr/2.
          endif

       CASE("S")

          if(F_RimeF<=5.0)then
             RHOX=100.
             if(NLICE>0.) then
                effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
             endif
          endif

       CASE("G")

          if(F_RimeF>5.0.and.F_RimeF<=20.0)then
             RHOX=400.
             if(NLICE>0.) then
                effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
             endif
          endif

       CASE("H")

          if(F_RimeF>20.0)then
             RHOX=900.
             if(NLICE>0.) then
                effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
             endif
          endif

       END SELECT


    endif

    !-----------------------------------------
    ! DIAMETER -> RADIUS
    !-----------------------------------------

    effr = 0.5*effr


  end function EFFR

!-------------------------------------------------------------------------------
!> @brief GAMMLN
!>
!> @param[in] XX
!> @return GAMMLN Returns the value of LN(GAMMA(XX)) FOR XX > 0.
!-------------------------------------------------------------------------------

  REAL FUNCTION GAMMLN(XX)
    !     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
    IMPLICIT NONE
    REAL, INTENT(IN):: XX
    DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
    DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
         COF = (/76.18009172947146D0, -86.50532032941677D0, &
         24.01409824083091D0, -1.231739572450155D0, &
         .1208650973866179D-2, -.5395239384953D-5/)
    DOUBLE PRECISION:: SER,TMP,X,Y
    INTEGER:: J

    X=XX
    Y=X
    TMP=X+5.5D0
    TMP=(X+0.5D0)*LOG(TMP)-TMP
    SER=1.000000000190015D0
    DO J=1,6
       Y=Y+1.D0
       SER=SER+COF(J)/Y
    ENDDO
    GAMMLN=TMP+LOG(STP*SER/X)
  END FUNCTION GAMMLN

!-------------------------------------------------------------------------------
!> @brief WGAMMA
!>
!> @param[in] y
!> @return WGAMMA
!-------------------------------------------------------------------------------

  REAL FUNCTION WGAMMA(y)

    IMPLICIT NONE
    REAL, INTENT(IN):: y

    ! real    :: GAMMLN

    WGAMMA = EXP(GAMMLN(y))

  END FUNCTION WGAMMA

  
end module mpas_crtm_diagnostics
