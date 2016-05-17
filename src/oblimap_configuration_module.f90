! File name: oblimap_configuration_module.f90          
! 
! Copyright (C) 2009 Thomas Reerink & Michael Kliphuis. This program
! is distributed under the terms of the GNU General Public License.
!
! This file is part of OBLIMAP.
!
! OBLIMAP is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! OBLIMAP is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License    
! along with OBLIMAP.  If not, see <http://www.gnu.org/licenses/>.
! 
!
! OBLIMAP is maintained by:
!
! Thomas Reerink
! Institute for Marine and Atmospheric Research Utrecht (IMAU)
! Utrecht University
! Princetonplein 5
! 3584 CC Utrecht
! The Netherlands
!
! email: <t.reerink@uu.nl>
!
! OBLIMAP is hosted on NeSCForge: NOT YET
!
! http://forge.nesc.ac.uk/projects/oblimap/
!

MODULE oblimap_configuration_module
      IMPLICIT NONE

      ! PRECISION
      ! =========
      ! The kind of real numbers used by default throughout the program.
      ! Reals should be declared as:
      !
      ! REAL(dp) :: example
      !  or
      ! REAL(KIND=dp) :: example
      !
      ! dp must be a PARAMETER
      INTEGER, PARAMETER :: dp  = KIND(1.0D0)  ! Kind of double precision numbers.

      
   ! CONFIG VARIABLES:
   !==================
   ! Variablles which are eventually set by the read_config_file subroutine:
      
      ! GRID SIZE AND SPACING
      ! =====================
      ! Number of grid points NX and NY and the spacing dx and dy (in m) in x and y-direction, number of grid points for the 
      ! linear vertical axis NZL (in the positive direction) and the stepsize z_step (in m) for making a cross section:
      INTEGER  :: NX_config     =   281                                      ! CONFIG variable
      INTEGER  :: NY_config     =   281                                      ! CONFIG variable  
      REAL(dp) :: dx_config     = 20000.0_dp                                 ! CONFIG variable
      REAL(dp) :: dy_config     = 20000.0_dp                                 ! CONFIG variable
      INTEGER  :: NZL_config    =   251                                      ! CONFIG variable
      INTEGER  :: z_step_config =    40                                      ! CONFIG variable

      ! Number of grid points in vertical direction:
      INTEGER  :: NZ_config     =    15                                      ! CONFIG variable
      ! Relative grid spacing in vertical direction. If k is a counter through the vertical layers, k=1 
      ! at the surface corresponding with zeta=0, and k=NZ at the bottom of the ice sheet corresponding with zeta=1:
      ! This CONFIG variable zeta_config is declared as a large array, because fortran does not allow a CONFIG/NAMELIST
      ! variable which is ALLOCATABLE, only the C%NZ first elements of this array will be used and have to be specified
      ! in the CONFIG file. 
      REAL(dp), DIMENSION(210), SAVE :: zeta_config = &
       (/0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
         0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp  /)
                                            
      
      ! GCM grid sizes
      ! ==============
      ! Number of GCM grid points in x and y-direction:
      INTEGER             :: NLON_config = 128                                                     ! CONFIG variable
      INTEGER             :: NLAT_config = 64                                                      ! CONFIG variable
                                                                                                   
      ! Number of extended grid points in x and y-direction of IM grid:                            
      INTEGER             :: NX_ext_config = 10                                                    ! CONFIG variable
      INTEGER             :: NY_ext_config = 10                                                    ! CONFIG variable

      ! MAPPING
      ! =======
      REAL(dp)           :: oblimap_allocate_factor_config    = 2                                  ! CONFIG variable
      CHARACTER(LEN=256) :: choice_projection_method_config   = 'oblique_stereographic_projection' ! CONFIG variable
      REAL(dp)           :: lambda_M_config                   = 320._dp                            ! CONFIG variable
      REAL(dp)           :: phi_M_config                      = 72._dp                             ! CONFIG variable
      REAL(dp)           :: alpha_stereographic_config        = 7.5_dp                             ! CONFIG variable
      INTEGER            :: shephard_exponent_quadrant_config = 2                                  ! CONFIG variable
      INTEGER            :: shephard_exponent_radius_config   = 2                                  ! CONFIG variable
      CHARACTER(LEN=256) :: gcm_input_filename_config         = &                                  ! CONFIG variable
                            './data/ccsm/ccsm_example_dec_feb_pd.nc'
      INTEGER            :: starting_recordnr_config = 1                                           ! CONFIG variable
      LOGICAL            :: choice_1D_lonlat_in_netcdf_config = .TRUE.                             ! CONFIG variable
      CHARACTER(LEN=128) :: gcm_x_axis_name_config            = 'lon'                              ! CONFIG variable
      CHARACTER(LEN=128) :: gcm_y_axis_name_config            = 'lat'                              ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_gcm_field_1_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_2_config         = 'Accumulation'                     ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_3_config         = 'PHIS'                             ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_4_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_5_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_6_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_7_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_8_config         = 'TS'                               ! CONFIG variable                
      CHARACTER(LEN=128) :: mapped_gcm_field_9_config         = 'TS'                               ! CONFIG variable                
      REAL(dp)           :: gcm_to_im_factor_field_1_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_2_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_3_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_4_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_5_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_6_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_7_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_8_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_factor_field_9_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_1_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_2_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_3_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_4_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_5_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_6_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_7_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_8_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: gcm_to_im_shift_field_9_config    = 0._dp                              ! CONFIG variable
      CHARACTER(LEN=256) :: im_created_filename_config        = 'created_IM.nc'                    ! CONFIG variable
      LOGICAL            :: choice_fast_map_GGM_to_IM_config  = .FALSE.                            ! CONFIG variable
      CHARACTER(LEN=256) :: input_fast_map_GCM_to_IM_config   = 'input_fast_GCM_to_IM_mapping.txt' ! CONFIG variable
      LOGICAL            :: choice_radius_GCM_to_IM_config    = .FALSE.                            ! CONFIG variable
      INTEGER            :: R_search_GCM_to_IM_config         = 16000                              ! CONFIG variable
      CHARACTER(LEN=256) :: im_input_filename_config          = 'created_IM.nc'                    ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_1_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_2_config          = 'MB_surface'                       ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_3_config          = 'Hs'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_4_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_5_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_6_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_7_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_8_config          = 'Ts'                               ! CONFIG variable
      CHARACTER(LEN=128) :: mapped_im_field_9_config          = 'Ts'                               ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_1_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_2_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_3_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_4_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_5_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_6_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_7_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_8_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_factor_field_9_config   = 1._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_1_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_2_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_3_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_4_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_5_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_6_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_7_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_8_config    = 0._dp                              ! CONFIG variable
      REAL(dp)           :: im_to_gcm_shift_field_9_config    = 0._dp                              ! CONFIG variable
      CHARACTER(LEN=256) :: gcm_created_filename_config       = 'to_and_fro_mapped_GCM.nc'         ! CONFIG variable
      LOGICAL            :: choice_fast_map_IM_to_GCM_config  = .FALSE.                            ! CONFIG variable
      CHARACTER(LEN=256) :: input_fast_map_IM_to_GCM_config   = 'input_fast_IM_to_GCM_mapping.txt' ! CONFIG variable
      LOGICAL            :: choice_radius_IM_to_GCM_config    = .TRUE.                             ! CONFIG variable
      INTEGER            :: R_search_IM_to_GCM_config         = 125000                             ! CONFIG variable
      LOGICAL            :: do_oblimap_post_processing_config = .FALSE.                            ! CONFIG variable

      ! The Earth Radius [meter]
      REAL(dp)           :: R_earth_config                    = 6.371221E6_dp                      ! CONFIG variable

    ! TYPE DEFENITION
    !================

      ! This TYPE contains all the information once the CONFIG file is read never will change during the run of the program 
      TYPE constants_type                
                                          ! GRID SIZES AND GRID SPACING
                                          !============================
        INTEGER                             :: NX          ! Number of grid points in x-direction
        INTEGER                             :: NY          ! Number of grid points in y-direction
        INTEGER                             :: NZ          ! Number of grid points in vertical direction for thermodynamics in ice sheet
        INTEGER                             :: NZL         ! Number of grid points in the Linear vertical direction
        REAL(dp)                            :: dx          ! Grid spacing in x-direction      
        REAL(dp)                            :: dy          ! Grid spacing in y-direction       
        REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta        ! Grid spacing in zeta-direction (zeta is the vertical scaled coordinate)  
        REAL(dp)                            :: z_step      ! Linear vertical grid size     
        INTEGER                             :: NLON        ! NX for GCM (longitudinal direction)
        INTEGER                             :: NLAT        ! NY for GCM (lateral direction) 
        INTEGER                             :: NX_ext      ! additional extended margin in x-direction
        INTEGER                             :: NY_ext      ! additional extended margin in y-direction
        INTEGER                             :: NX_EXTENDED ! Total size of the extended grid in x-direction
        INTEGER                             :: NY_EXTENDED ! Total size of the extended grid in y-direction

                                          ! MATHEMATICAL CONSTANTS
                                          !=======================
        REAL(dp)                            :: pi
        REAL(dp)                            :: deg2rad  ! Conversion factor between radians and degrees
        REAL(dp)                            :: rad2deg  ! Conversion factor between degrees and radians

                                          ! MAPPING
                                          !========
        REAL(dp)                            :: oblimap_allocate_factor
        CHARACTER(LEN=256)                  :: choice_projection_method
        REAL(dp)                            :: lambda_M
        REAL(dp)                            :: phi_M
        REAL(dp)                            :: alpha_stereographic
        INTEGER                             :: shephard_exponent_quadrant
        INTEGER                             :: shephard_exponent_radius
        CHARACTER(LEN=128)                  :: gcm_input_filename
        INTEGER                             :: starting_recordnr
        LOGICAL                             :: choice_1D_lonlat_in_netcdf
        CHARACTER(LEN=128)                  :: gcm_x_axis_name
        CHARACTER(LEN=128)                  :: gcm_y_axis_name
        CHARACTER(LEN=128)                  :: mapped_gcm_field_1
        CHARACTER(LEN=128)                  :: mapped_gcm_field_2
        CHARACTER(LEN=128)                  :: mapped_gcm_field_3
        CHARACTER(LEN=128)                  :: mapped_gcm_field_4
        CHARACTER(LEN=128)                  :: mapped_gcm_field_5
        CHARACTER(LEN=128)                  :: mapped_gcm_field_6
        CHARACTER(LEN=128)                  :: mapped_gcm_field_7
        CHARACTER(LEN=128)                  :: mapped_gcm_field_8
        CHARACTER(LEN=128)                  :: mapped_gcm_field_9
        REAL(dp)                            :: gcm_to_im_factor_field_1
        REAL(dp)                            :: gcm_to_im_factor_field_2
        REAL(dp)                            :: gcm_to_im_factor_field_3
        REAL(dp)                            :: gcm_to_im_factor_field_4
        REAL(dp)                            :: gcm_to_im_factor_field_5
        REAL(dp)                            :: gcm_to_im_factor_field_6
        REAL(dp)                            :: gcm_to_im_factor_field_7
        REAL(dp)                            :: gcm_to_im_factor_field_8
        REAL(dp)                            :: gcm_to_im_factor_field_9
        REAL(dp)                            :: gcm_to_im_shift_field_1
        REAL(dp)                            :: gcm_to_im_shift_field_2
        REAL(dp)                            :: gcm_to_im_shift_field_3
        REAL(dp)                            :: gcm_to_im_shift_field_4
        REAL(dp)                            :: gcm_to_im_shift_field_5
        REAL(dp)                            :: gcm_to_im_shift_field_6
        REAL(dp)                            :: gcm_to_im_shift_field_7
        REAL(dp)                            :: gcm_to_im_shift_field_8
        REAL(dp)                            :: gcm_to_im_shift_field_9
        CHARACTER(LEN=128)                  :: im_created_filename
        LOGICAL                             :: choice_fast_map_GGM_to_IM
        CHARACTER(LEN=128)                  :: input_fast_map_GCM_to_IM
        LOGICAL                             :: choice_radius_GCM_to_IM
        INTEGER                             :: R_search_GCM_to_IM
        CHARACTER(LEN=128)                  :: im_input_filename
        CHARACTER(LEN=128)                  :: mapped_im_field_1
        CHARACTER(LEN=128)                  :: mapped_im_field_2
        CHARACTER(LEN=128)                  :: mapped_im_field_3
        CHARACTER(LEN=128)                  :: mapped_im_field_4
        CHARACTER(LEN=128)                  :: mapped_im_field_5
        CHARACTER(LEN=128)                  :: mapped_im_field_6
        CHARACTER(LEN=128)                  :: mapped_im_field_7
        CHARACTER(LEN=128)                  :: mapped_im_field_8
        CHARACTER(LEN=128)                  :: mapped_im_field_9
        REAL(dp)                            :: im_to_gcm_factor_field_1
        REAL(dp)                            :: im_to_gcm_factor_field_2
        REAL(dp)                            :: im_to_gcm_factor_field_3
        REAL(dp)                            :: im_to_gcm_factor_field_4
        REAL(dp)                            :: im_to_gcm_factor_field_5
        REAL(dp)                            :: im_to_gcm_factor_field_6
        REAL(dp)                            :: im_to_gcm_factor_field_7
        REAL(dp)                            :: im_to_gcm_factor_field_8
        REAL(dp)                            :: im_to_gcm_factor_field_9
        REAL(dp)                            :: im_to_gcm_shift_field_1
        REAL(dp)                            :: im_to_gcm_shift_field_2
        REAL(dp)                            :: im_to_gcm_shift_field_3
        REAL(dp)                            :: im_to_gcm_shift_field_4
        REAL(dp)                            :: im_to_gcm_shift_field_5
        REAL(dp)                            :: im_to_gcm_shift_field_6
        REAL(dp)                            :: im_to_gcm_shift_field_7
        REAL(dp)                            :: im_to_gcm_shift_field_8
        REAL(dp)                            :: im_to_gcm_shift_field_9
        CHARACTER(LEN=128)                  :: gcm_created_filename
        LOGICAL                             :: choice_fast_map_IM_to_GCM
        CHARACTER(LEN=128)                  :: input_fast_map_IM_to_GCM
        LOGICAL                             :: choice_radius_IM_to_GCM
        INTEGER                             :: R_search_IM_to_GCM
        LOGICAL                             :: do_oblimap_post_processing

        REAL(dp)                            :: R_earth
        
                                          ! MAPPING ELLIPSOID
                                          !==================
        REAL(dp)                            :: a
        REAL(dp)                            :: e
        REAL(dp)                            :: am
        REAL(dp)                            :: akm
        REAL(dp)                            :: chi_M
        
        REAL(dp)                            :: q_M      
        REAL(dp)                            :: q_polar  
        REAL(dp)                            :: beta_M   
        REAL(dp)                            :: R_q_polar
        REAL(dp)                            :: D        
      END TYPE constants_type
      
      ! C is the 'struct' containing all the Constants from the CONFIG file and/or the defaults
      TYPE(constants_type), SAVE :: C



CONTAINS
  SUBROUTINE read_config_file()
    ! This subroutine reads the CONFIG variables from a configuration file. The name of the
    ! configuration file should be specified on the command line. If no name is specified on
    ! the command line, then the default values, as specified in this module, are used.
    IMPLICIT NONE
    
    ! Local variables:
    CHARACTER(LEN=128) :: config_filename
    INTEGER, PARAMETER :: config_unit = 28 ! Unit number which is used for the configuration file.
    INTEGER            :: ios
    INTEGER            :: iargc

     ! List of items in the configuration file:
     NAMELIST /CONFIG/NLON_config                               , &
                      NLAT_config                               , &
                      NX_config                                 , &
                      NY_config                                 , &
                      dx_config                                 , &
                      dy_config                                 , &
                      NX_ext_config                             , &
                      NY_ext_config                             , &
                      oblimap_allocate_factor_config            , &
                      choice_projection_method_config           , &
                      lambda_M_config                           , &
                      phi_M_config                              , &
                      alpha_stereographic_config                , &
                      shephard_exponent_quadrant_config         , &
                      shephard_exponent_radius_config           , &
                      gcm_input_filename_config                 , &
                      starting_recordnr_config                  , &
                      choice_1D_lonlat_in_netcdf_config         , &
                      gcm_x_axis_name_config                    , &
                      gcm_y_axis_name_config                    , &
                      mapped_gcm_field_1_config                 , &
                      mapped_gcm_field_2_config                 , &
                      mapped_gcm_field_3_config                 , &
                      mapped_gcm_field_4_config                 , &
                      mapped_gcm_field_5_config                 , &
                      mapped_gcm_field_6_config                 , &
                      mapped_gcm_field_7_config                 , &
                      mapped_gcm_field_8_config                 , &
                      mapped_gcm_field_9_config                 , &
                      gcm_to_im_factor_field_1_config           , &
                      gcm_to_im_factor_field_2_config           , &
                      gcm_to_im_factor_field_3_config           , &
                      gcm_to_im_factor_field_4_config           , &
                      gcm_to_im_factor_field_5_config           , &
                      gcm_to_im_factor_field_6_config           , &
                      gcm_to_im_factor_field_7_config           , &
                      gcm_to_im_factor_field_8_config           , &
                      gcm_to_im_factor_field_9_config           , &
                      gcm_to_im_shift_field_1_config            , &
                      gcm_to_im_shift_field_2_config            , &
                      gcm_to_im_shift_field_3_config            , &
                      gcm_to_im_shift_field_4_config            , &
                      gcm_to_im_shift_field_5_config            , &
                      gcm_to_im_shift_field_6_config            , &
                      gcm_to_im_shift_field_7_config            , &
                      gcm_to_im_shift_field_8_config            , &
                      gcm_to_im_shift_field_9_config            , &
                      im_created_filename_config                , &
                      choice_fast_map_GGM_to_IM_config          , &
                      input_fast_map_GCM_to_IM_config           , &
                      choice_radius_GCM_to_IM_config            , &
                      R_search_GCM_to_IM_config                 , &
                      im_input_filename_config                  , &
                      mapped_im_field_1_config                  , &
                      mapped_im_field_2_config                  , &
                      mapped_im_field_3_config                  , &
                      mapped_im_field_4_config                  , &
                      mapped_im_field_5_config                  , &
                      mapped_im_field_6_config                  , &
                      mapped_im_field_7_config                  , &
                      mapped_im_field_8_config                  , &
                      mapped_im_field_9_config                  , &
                      im_to_gcm_factor_field_1_config           , &
                      im_to_gcm_factor_field_2_config           , &
                      im_to_gcm_factor_field_3_config           , &
                      im_to_gcm_factor_field_4_config           , &
                      im_to_gcm_factor_field_5_config           , &
                      im_to_gcm_factor_field_6_config           , &
                      im_to_gcm_factor_field_7_config           , &
                      im_to_gcm_factor_field_8_config           , &
                      im_to_gcm_factor_field_9_config           , &
                      im_to_gcm_shift_field_1_config            , &
                      im_to_gcm_shift_field_2_config            , &
                      im_to_gcm_shift_field_3_config            , &
                      im_to_gcm_shift_field_4_config            , &
                      im_to_gcm_shift_field_5_config            , &
                      im_to_gcm_shift_field_6_config            , &
                      im_to_gcm_shift_field_7_config            , &
                      im_to_gcm_shift_field_8_config            , &
                      im_to_gcm_shift_field_9_config            , &
                      gcm_created_filename_config               , &
                      choice_fast_map_IM_to_GCM_config          , &
                      input_fast_map_IM_to_GCM_config           , &
                      choice_radius_IM_to_GCM_config            , &
                      R_search_IM_to_GCM_config                 , &
                      do_oblimap_post_processing_config         , &
                      R_earth_config
         
    SELECT CASE(iargc())  
    CASE(0)
     ! Don't read the configuration file in this case, just use the default values.
    CASE(1)
     ! Get the name of the configuration file, open this file and read it:
     CALL getarg(1, config_filename)
     OPEN(UNIT=config_unit, FILE=TRIM(config_filename), STATUS='OLD', ACTION='READ', iostat=ios)
     IF(ios /= 0) THEN
      WRITE(UNIT=*, FMT=*) 'Could not open configuratio file: ', TRIM(config_filename)
      STOP
     END IF
     ! In the following statement the entire configuration file is read, using the namelist (NML=CONFIG)
     READ(UNIT=config_unit, NML=CONFIG,IOSTAT=ios)
     IF(ios /= 0) THEN
      WRITE(UNIT=*, FMT=*) 'Error while reading configuration file: ', TRIM(config_filename)
      STOP
     END IF
     CLOSE(UNIT=config_unit)
    CASE DEFAULT
     ! Specifying more than one command line argument is not allowed. Print an error message:
     WRITE(UNIT=*, FMT='(A)') 'Error: Too many command line arguments!'
     STOP
    END SELECT
  END SUBROUTINE read_config_file



  SUBROUTINE initialize_constants()
    ! This routine puts all the constants which will never change during the run after the CONFIG file
    ! has been read, into a special constant 'struct' 
    IMPLICIT NONE

    ! Number of grid points in x and y-direction and grid spacing in meters for x and y-direction:
    C%NX     = NX_config
    C%dx     = dx_config
    C%NY     = NY_config
    C%dy     = dy_config
 
    ! Number of grid points in vertical direction:
    C%NZ     = NZ_config
    ALLOCATE(C%zeta(C%NZ))                 
    C%zeta   = zeta_config(1:C%NZ) ! Fortran does not allow a CONFIG/NAMELIST variable to be ALLOCATABLE, therefore this way
  
    ! Number of grid points in the linear vertical direction and the grid spacing:
    C%NZL    = NZL_config
    C%z_step = z_step_config 

    ! GCM grid sizes:
    C%NLON   = NLON_config
    C%NLAT   = NLAT_config
    C%NX_ext = NX_ext_config
    C%NY_ext = NY_ext_config
    C%NX_EXTENDED = C%NX + 2 * C%NX_ext
    C%NY_EXTENDED = C%NY + 2 * C%NY_ext

    C%pi      = 2._dp*ACOS(0._dp)      ! Just pi=3.14159... exactly
    C%deg2rad = C%pi/180._dp           ! Conversion factor between radians and degrees
    C%rad2deg = 180._dp/C%pi           ! Conversion factor between degrees and radians

    C%oblimap_allocate_factor  = oblimap_allocate_factor_config
    C%choice_projection_method = choice_projection_method_config
    
    ! Assign a C%lambda_M value for North and South Pole which generate the correct IM coordinate system, 
    ! see Reerink et al. (2010) equation (2.3) or equation (A.53):
    IF(phi_M_config == -90._dp .OR. phi_M_config == 90._dp) lambda_M_config = 0._dp
    ! Coordinate (C%lamda_M, C%phi_M) is the middle of the GCM's longitude-latitude region of interest that 
    ! will be projected to the IM, convert the degrees to radians:
    C%lambda_M = C%deg2rad * lambda_M_config  
    C%phi_M    = C%deg2rad * phi_M_config
        
    ! The exact projection plane can be adjusted by specifying a certain angle alpha_stereographic
    ! This projection plane below and parallel to the tangent plane in (C%lamda_M, C%phi_M). 
    C%alpha_stereographic = C%deg2rad * alpha_stereographic_config 

    C%shephard_exponent_quadrant = shephard_exponent_quadrant_config
    C%shephard_exponent_radius   = shephard_exponent_radius_config
    C%gcm_input_filename         = gcm_input_filename_config
    C%starting_recordnr          = starting_recordnr_config 
    C%choice_1D_lonlat_in_netcdf = choice_1D_lonlat_in_netcdf_config
    C%gcm_x_axis_name            = gcm_x_axis_name_config
    C%gcm_y_axis_name            = gcm_y_axis_name_config
    C%mapped_gcm_field_1         = mapped_gcm_field_1_config
    C%mapped_gcm_field_2         = mapped_gcm_field_2_config
    C%mapped_gcm_field_3         = mapped_gcm_field_3_config
    C%mapped_gcm_field_4         = mapped_gcm_field_4_config
    C%mapped_gcm_field_5         = mapped_gcm_field_5_config
    C%mapped_gcm_field_6         = mapped_gcm_field_6_config
    C%mapped_gcm_field_7         = mapped_gcm_field_7_config
    C%mapped_gcm_field_8         = mapped_gcm_field_8_config
    C%mapped_gcm_field_9         = mapped_gcm_field_9_config
    C%gcm_to_im_factor_field_1   = gcm_to_im_factor_field_1_config   ! The factor from GCM to IM concerning field_1
    C%gcm_to_im_factor_field_2   = gcm_to_im_factor_field_2_config   ! The factor from GCM to IM concerning field_2
    C%gcm_to_im_factor_field_3   = gcm_to_im_factor_field_3_config   ! The factor from GCM to IM concerning field_3
    C%gcm_to_im_factor_field_4   = gcm_to_im_factor_field_4_config   ! The factor from GCM to IM concerning field_4
    C%gcm_to_im_factor_field_5   = gcm_to_im_factor_field_5_config   ! The factor from GCM to IM concerning field_5
    C%gcm_to_im_factor_field_6   = gcm_to_im_factor_field_6_config   ! The factor from GCM to IM concerning field_6
    C%gcm_to_im_factor_field_7   = gcm_to_im_factor_field_7_config   ! The factor from GCM to IM concerning field_7
    C%gcm_to_im_factor_field_8   = gcm_to_im_factor_field_8_config   ! The factor from GCM to IM concerning field_8
    C%gcm_to_im_factor_field_9   = gcm_to_im_factor_field_9_config   ! The factor from GCM to IM concerning field_9
    C%gcm_to_im_shift_field_1    = gcm_to_im_shift_field_1_config    ! The shift  from GCM to IM concerning field_1 
    C%gcm_to_im_shift_field_2    = gcm_to_im_shift_field_2_config    ! The shift  from GCM to IM concerning field_2 
    C%gcm_to_im_shift_field_3    = gcm_to_im_shift_field_3_config    ! The shift  from GCM to IM concerning field_3 
    C%gcm_to_im_shift_field_4    = gcm_to_im_shift_field_4_config    ! The shift  from GCM to IM concerning field_4 
    C%gcm_to_im_shift_field_5    = gcm_to_im_shift_field_5_config    ! The shift  from GCM to IM concerning field_5 
    C%gcm_to_im_shift_field_6    = gcm_to_im_shift_field_6_config    ! The shift  from GCM to IM concerning field_6 
    C%gcm_to_im_shift_field_7    = gcm_to_im_shift_field_7_config    ! The shift  from GCM to IM concerning field_7 
    C%gcm_to_im_shift_field_8    = gcm_to_im_shift_field_8_config    ! The shift  from GCM to IM concerning field_8 
    C%gcm_to_im_shift_field_9    = gcm_to_im_shift_field_9_config    ! The shift  from GCM to IM concerning field_9 
    C%im_created_filename        = im_created_filename_config 
    C%choice_fast_map_GGM_to_IM  = choice_fast_map_GGM_to_IM_config
    C%input_fast_map_GCM_to_IM   = input_fast_map_GCM_to_IM_config
    C%choice_radius_GCM_to_IM    = choice_radius_GCM_to_IM_config
    C%R_search_IM_to_GCM         = R_search_IM_to_GCM_config
    C%im_input_filename          = im_input_filename_config
    C%mapped_im_field_1          = mapped_im_field_1_config
    C%mapped_im_field_2          = mapped_im_field_2_config
    C%mapped_im_field_3          = mapped_im_field_3_config
    C%mapped_im_field_4          = mapped_im_field_4_config
    C%mapped_im_field_5          = mapped_im_field_5_config
    C%mapped_im_field_6          = mapped_im_field_6_config
    C%mapped_im_field_7          = mapped_im_field_7_config
    C%mapped_im_field_8          = mapped_im_field_8_config
    C%mapped_im_field_9          = mapped_im_field_9_config
    C%im_to_gcm_factor_field_1   = im_to_gcm_factor_field_1_config   ! The factor from IM to GCM concerning field_1
    C%im_to_gcm_factor_field_2   = im_to_gcm_factor_field_2_config   ! The factor from IM to GCM concerning field_2
    C%im_to_gcm_factor_field_3   = im_to_gcm_factor_field_3_config   ! The factor from IM to GCM concerning field_3
    C%im_to_gcm_factor_field_4   = im_to_gcm_factor_field_4_config   ! The factor from IM to GCM concerning field_4
    C%im_to_gcm_factor_field_5   = im_to_gcm_factor_field_5_config   ! The factor from IM to GCM concerning field_5
    C%im_to_gcm_factor_field_6   = im_to_gcm_factor_field_6_config   ! The factor from IM to GCM concerning field_6
    C%im_to_gcm_factor_field_7   = im_to_gcm_factor_field_7_config   ! The factor from IM to GCM concerning field_7
    C%im_to_gcm_factor_field_8   = im_to_gcm_factor_field_8_config   ! The factor from IM to GCM concerning field_8
    C%im_to_gcm_factor_field_9   = im_to_gcm_factor_field_9_config   ! The factor from IM to GCM concerning field_9
    C%im_to_gcm_shift_field_1    = im_to_gcm_shift_field_1_config    ! The shift  from IM to GCM concerning field_1 
    C%im_to_gcm_shift_field_2    = im_to_gcm_shift_field_2_config    ! The shift  from IM to GCM concerning field_2 
    C%im_to_gcm_shift_field_3    = im_to_gcm_shift_field_3_config    ! The shift  from IM to GCM concerning field_3 
    C%im_to_gcm_shift_field_4    = im_to_gcm_shift_field_4_config    ! The shift  from IM to GCM concerning field_4 
    C%im_to_gcm_shift_field_5    = im_to_gcm_shift_field_5_config    ! The shift  from IM to GCM concerning field_5 
    C%im_to_gcm_shift_field_6    = im_to_gcm_shift_field_6_config    ! The shift  from IM to GCM concerning field_6 
    C%im_to_gcm_shift_field_7    = im_to_gcm_shift_field_7_config    ! The shift  from IM to GCM concerning field_7 
    C%im_to_gcm_shift_field_8    = im_to_gcm_shift_field_8_config    ! The shift  from IM to GCM concerning field_8 
    C%im_to_gcm_shift_field_9    = im_to_gcm_shift_field_9_config    ! The shift  from IM to GCM concerning field_9 
    C%gcm_created_filename       = gcm_created_filename_config
    C%choice_fast_map_IM_to_GCM  = choice_fast_map_IM_to_GCM_config
    C%input_fast_map_IM_to_GCM   = input_fast_map_IM_to_GCM_config
    C%choice_radius_IM_to_GCM    = choice_radius_IM_to_GCM_config
    C%R_search_GCM_to_IM         = R_search_GCM_to_IM_config
    C%do_oblimap_post_processing = do_oblimap_post_processing_config

    C%R_earth                    = R_earth_config

    C%a      = 6378137._dp    ! equatorial ellipsoid radius, a in Snyder, WGS84
    C%e      = 0.081819191_dp ! excentricity of the ellipsoid, e in Snyder, WGS84, see Snyder p. 13
    ! See equations (14-15) and (21-27) on page 160 in Snyder (1987), akm corresponds with 2a*k0*m_1 in Snyder:
    C%am     = C%a * (COS(C%phi_M) / SQRT(1._dp - (C%e * SIN(C%phi_M))**2))
    C%akm    = (1._dp + COS(C%alpha_stereographic)) * C%am
    ! See equations (3-1a) on page 160 in Snyder (1987),  chi_M corresponds with chi_1 in Snyder:
    C%chi_M  = 2._dp * ATAN(SQRT(((1._dp +       SIN(C%phi_M)) / (1._dp -       SIN(C%phi_M))) * &
                                 ((1._dp - C%e * SIN(C%phi_M)) / (1._dp + C%e * SIN(C%phi_M)))**(C%e))) - 0.5_dp * C%pi

    ! See equation (3-12) on page 187 in Snyder (1987):
    C%q_M       = (1._dp - C%e**2) * ((SIN(C%phi_M) / (1._dp - (C%e * SIN(C%phi_M))**2)) - (1._dp / (2._dp * C%e)) * LOG((1._dp - C%e * SIN(C%phi_M)) / (1._dp + C%e * SIN(C%phi_M)))) 
    ! See equation (3-12) on page 187 in Snyder (1987):
    C%q_polar   = (1._dp - C%e**2) * ((1._dp / (1._dp - C%e**2)) - (1._dp / (2._dp * C%e)) * LOG((1._dp - C%e) / (1._dp + C%e)))
    ! See equation (3-11) on page 187 in Snyder (1987):
    C%beta_M    = ASIN(C%q_M / C%q_polar)
    ! See equation (3-13) on page 187 in Snyder (1987):
    C%R_q_polar = C%a * SQRT(0.5_dp * C%q_polar)
    ! See equation (24-20) on page 187 in Snyder (1987):
    C%D         = C%am / (C%R_q_polar * COS(C%phi_M))
  END SUBROUTINE initialize_constants  



  SUBROUTINE finalize_constants()
    DEALLOCATE(C%zeta)                 
  END SUBROUTINE finalize_constants

END MODULE oblimap_configuration_module
