! File name: oblimap_im_to_gcm_program.f90          
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

PROGRAM oblimap_im_to_gcm_program
  ! This program reads several field variables from a rectangular Ice Model (IM) netcdf file which 
  ! are projected by an inverse oblique stereographic projection on a Global Circulation Model (GCM)
  ! grid. The IM grid point coordinates are projected on the GCM grid, where in general they fall
  ! irregular and between the GCM grid points. Therefore the field values defined on the IM grid  
  ! points have to be interpolated at each GCM grid point, with help of the nearby projected IM
  ! grid points. Two interpolation methods are available in this program both based on a distance 
  ! weigthing Shepard technique. One method, the 'quadrant method', searches within each quadrant 
  ! around each GCM grid point the nearest projected IM point and interpolates with help of this
  ! four points to obtain the field value at such a GCM grid point. Another method, the 'radius
  ! method', searches all projected IM points within a certain radius and interpolates with help 
  ! of these points to obtain the field value at such a GCM grid point. Finally the projected and
  ! interpolated GCM fields are written to an GCM netcdf file. The GCM points which were not affected
  ! by the mapping (= projection + interpolation) remain their values as before.
  !    
  ! The GCM geographical coordinates are in longitude (lon) and latitude (lat) in degrees, while the 
  ! IM rectangular coordinates are in x and y in meters. The (lon,lat) coordinates are defined in the
  ! curved spherical surface S, and the (x,y) coordinates in the flat surface S'. For a more extended
  ! description of the projection and the interpolation method see:
  !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
  !    
  ! In the NAMELIST or CONFIG file, several options can be changed without compiling the program.
  ! For example the center of the area of interest (the center of the projected area) can be specified
  ! by setting the:
  !   lamda_M_config = 320  (values for the case of Greenland)
  !   phi_M_config   = 72   (values for the case of Greenland)
  ! in the CONFIG file. The exact (inverse) oblique stereographic projection can be chosen with
  !   alpha_stereographic_config = 19 (a not too bad example)
  ! The center of the IM grid will coincide with (lamda_M_config,phi_M_config), and the extensions of
  ! the IM grid are determined by the IM grid spacings C%dx and C%dy and the IM grid sizes C%NX and C%NY.
  !
  !
  ! For the mapping from IM to GCM we make use of initial_gcm_field_*, created_gcm_field_*, and gcm_field_*.
  ! 
  ! initial_gcm_field_*
  !   Contains the initial GCM fields. This fields are used because a large part of each GCM field is 
  !   unaffected by the to and fro mapping because the IM grid is only local. The unaffected (mask(i,j) == 0)
  !   field values of the initial fields are afterwards merged with the to and fro mapped values. Besides
  !   in the 'scan' phase the longitude and latitude coordinates of the GCM grid points have to be read
  !   from the GCM file which contains these initial_gcm_field_* fields.
  ! 
  ! created_gcm_field_*
  !   Contains the GCM fields which are mapped (= projected + interpolated) from the IM grid on to the GCM
  !   grid, these fields only have non-zero values on the grid points (i,j) with mask(i,j) == 1.
  ! 
  ! gcm_field_*
  !   Contains the to and fro mapped GCM fields merged with the initial_gcm_field_*. 
  !
  USE oblimap_configuration_module, ONLY: dp, C, read_config_file, initialize_constants, finalize_constants
  USE oblimap_mapping_module, ONLY: fast_mapping
  USE oblimap_read_and_write_module, ONLY: read_data_gcm, read_data_gcm_with_2D_lonlat, read_data_im, &
    write_data_for_gcm, write_data_for_gcm_with_2D_lonlat 
  USE oblimap_scan_contributions_module, ONLY: make_mask, scan_with_radius_method_im_to_gcm, scan_with_quadrant_method_im_to_gcm
  USE oblimap_post_processing_module, ONLY: calculate_average_differences
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_1          ! A variable at the IM grid which will be projected to the GCM grid
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_2          ! A variable at the IM grid which will be projected to the GCM grid
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_3          ! A variable at the IM grid which will be projected to the GCM grid
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_1_extended   
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_2_extended   
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_3_extended   

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: initial_gcm_field_1    
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: initial_gcm_field_2   
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: initial_gcm_field_3   
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: created_gcm_field_1 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: created_gcm_field_2 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: created_gcm_field_3 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_1  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_2  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_3  
  INTEGER,  DIMENSION(:,:), ALLOCATABLE :: mask
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lon_gcm
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lat_gcm
  
  ! Read the config-variables which are declared in the grid module from the configuration file:
  CALL read_config_file()
  ! Initialization of the struckt  C% :
  CALL initialize_constants()

  ! Allocate all arrays for which the dimension sizes are given in the CONFIG file:
  ALLOCATE(im_field_1           (C%NX         ,C%NY         ))
  ALLOCATE(im_field_2           (C%NX         ,C%NY         ))
  ALLOCATE(im_field_3           (C%NX         ,C%NY         ))
  ALLOCATE(im_field_1_extended  (C%NX_EXTENDED,C%NY_EXTENDED))
  ALLOCATE(im_field_2_extended  (C%NX_EXTENDED,C%NY_EXTENDED))
  ALLOCATE(im_field_3_extended  (C%NX_EXTENDED,C%NY_EXTENDED))
  ALLOCATE(initial_gcm_field_1  (C%NLON       ,C%NLAT       ))
  ALLOCATE(initial_gcm_field_2  (C%NLON       ,C%NLAT       ))
  ALLOCATE(initial_gcm_field_3  (C%NLON       ,C%NLAT       ))
  ALLOCATE(created_gcm_field_1  (C%NLON       ,C%NLAT       ))
  ALLOCATE(created_gcm_field_2  (C%NLON       ,C%NLAT       ))
  ALLOCATE(created_gcm_field_3  (C%NLON       ,C%NLAT       ))
  ALLOCATE(gcm_field_1          (C%NLON       ,C%NLAT       ))
  ALLOCATE(gcm_field_2          (C%NLON       ,C%NLAT       ))
  ALLOCATE(gcm_field_3          (C%NLON       ,C%NLAT       ))
  ALLOCATE(mask                 (C%NLON       ,C%NLAT       ))
  ALLOCATE(lon_gcm              (C%NLON       ,C%NLAT       ))
  ALLOCATE(lat_gcm              (C%NLON       ,C%NLAT       ))

  ! Reading the GCM field variables 'C%mapped_gcm_field_*' and the longitude and latitude coordinates
  ! of the GCM grid from a GCM-netcdf file. The initial GCM fields are needed to merge the mapped results.
  IF(C%choice_1D_lonlat_in_netcdf .EQV. .TRUE.) THEN
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: lon_gcm, lat_gcm, initial_gcm_field_1, initial_gcm_field_2, initial_gcm_field_3
   CALL read_data_gcm(C%gcm_input_filename, C%mapped_gcm_field_1, lon_gcm, lat_gcm, initial_gcm_field_1, &
                                            C%mapped_gcm_field_2,                   initial_gcm_field_2, &
                                            C%mapped_gcm_field_3,                   initial_gcm_field_3)
  ELSE
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: lon_gcm, lat_gcm, initial_gcm_field_1, initial_gcm_field_2, initial_gcm_field_3
   CALL read_data_gcm_with_2D_lonlat(C%gcm_input_filename, C%mapped_gcm_field_1, lon_gcm, lat_gcm, initial_gcm_field_1, &
                                                           C%mapped_gcm_field_2,                   initial_gcm_field_2, &
                                                           C%mapped_gcm_field_3,                   initial_gcm_field_3)
  
   ! Adhoc racmo coordinate shift correction:
   WHERE(lon_gcm < 0._dp) lon_gcm = lon_gcm + 360._dp
  END IF

  ! From the IM netcdf file we read the IM variables C%mapped_im_field_1 -- C%mapped_im_field_3 which will be mapped. 
  ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
  ! Output: im_field_1, im_field_2, im_field_3
  CALL read_data_im(C%im_input_filename, C%mapped_im_field_1, im_field_1, &
                                         C%mapped_im_field_2, im_field_2, &
                                         C%mapped_im_field_3, im_field_3)
  im_field_1 = C%im_to_gcm_factor_field_1 * im_field_1 + C%im_to_gcm_shift_field_1
  im_field_2 = C%im_to_gcm_factor_field_2 * im_field_2 + C%im_to_gcm_shift_field_2
  im_field_3 = C%im_to_gcm_factor_field_3 * im_field_3 + C%im_to_gcm_shift_field_3
       
  ! Output: im_field_1_extended
  CALL initialize_extended_im_field(im_field_1, im_field_1_extended)
  ! Output: im_field_2_extended
  CALL initialize_extended_im_field(im_field_2, im_field_2_extended)
  ! Output: im_field_3_extended
  CALL initialize_extended_im_field(im_field_3, im_field_3_extended)

  ! Before the fast mapping routine can be used, the C%input_fast_map_IM_to_GCM file has to be created. This file
  ! contains for each target grid point the coordinates and the relative distances of the contributing points which 
  ! are used for the field interpolation. This contributing points of one such a target point are the nearest
  ! projected points selected by the radius method or by the quadrant method. This 'scan' part in fact contains the
  ! most technical and CPU consuming part of OBLIMAP covering the projection and the selection of points for the 
  ! interpolation (and keeping the relative distances of each projected point relative to the target point).
  IF(C%choice_fast_map_IM_to_GCM .EQV. .FALSE.) THEN
   IF(C%choice_radius_IM_to_GCM .EQV. .TRUE.) THEN
    ! Output: the C%input_fast_map_IM_to_GCM is created 
    CALL scan_with_radius_method_im_to_gcm(lon_gcm, lat_gcm)
   ELSE
    ! Output: the C%input_fast_map_IM_to_GCM is created 
    CALL scan_with_quadrant_method_im_to_gcm(lon_gcm, lat_gcm)
   END IF
  END IF

  ! The IM extended fields are mapped (= projected + interpolated) on the GCM grid. For each target grid point the coordinates and 
  ! the relative distances of the nearest projected points are stored in the C%input_fast_map_IM_to_GCM file, these are 
  ! used to map the fields with the fast mapping:
  ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
  ! Output: created_gcm_field_1, created_gcm_field_2, created_gcm_field_3
  CALL fast_mapping(C%input_fast_map_IM_to_GCM, C%NX_EXTENDED, C%NY_EXTENDED, C%NLON, C%NLAT, im_field_1_extended, created_gcm_field_1, &
                                                                                              im_field_2_extended, created_gcm_field_2, &
                                                                                              im_field_3_extended, created_gcm_field_3)
  
  ! The IM grid is defined on just a small part of the world (e.g. Antarctica, Greenland etc.) see the CONFIG file
  ! After to and fro mapping only field values for points within the IM domain are 'known'. 
  ! So all other points in the GCM field which did not participate in the mapping are taken equal 
  ! to the ones in the original pre-mapped GCM file:
  !  gcm_field(mask=0,t+1) = gcm_field(mask=0,t) 
  ! Output: mask
  CALL make_mask(lon_gcm, lat_gcm, mask) 
  
  ! Compose the fields which are given back to GCM. Part of these fields (the coordinates with mask = 0) are
  ! the same as the original GCM fields: initial_gcm_field_*, and the rest of them 
  ! equal the mapped values as in created_gcm_field_*. 
  gcm_field_1 = (1 - mask) * initial_gcm_field_1 + created_gcm_field_1
  gcm_field_2 = (1 - mask) * initial_gcm_field_2 + created_gcm_field_2
  gcm_field_3 = (1 - mask) * initial_gcm_field_3 + created_gcm_field_3

  IF(C%do_oblimap_post_processing .EQV. .TRUE.) THEN
   CALL calculate_average_differences(C%gcm_to_im_factor_field_1 * im_field_1 + C%gcm_to_im_shift_field_1, C%gcm_to_im_factor_field_1 * gcm_field_1 + C%gcm_to_im_shift_field_1, &
                                      C%gcm_to_im_factor_field_1 * initial_gcm_field_1 + C%gcm_to_im_shift_field_1, mask, C%mapped_im_field_1, lon_gcm, lat_gcm)
   CALL calculate_average_differences(C%gcm_to_im_factor_field_2 * im_field_2 + C%gcm_to_im_shift_field_2, C%gcm_to_im_factor_field_2 * gcm_field_2 + C%gcm_to_im_shift_field_2, &
                                      C%gcm_to_im_factor_field_2 * initial_gcm_field_2 + C%gcm_to_im_shift_field_2, mask, C%mapped_im_field_2, lon_gcm, lat_gcm)
   CALL calculate_average_differences(C%gcm_to_im_factor_field_3 * im_field_3 + C%gcm_to_im_shift_field_3, C%gcm_to_im_factor_field_3 * gcm_field_3 + C%gcm_to_im_shift_field_3, &
                                      C%gcm_to_im_factor_field_3 * initial_gcm_field_3 + C%gcm_to_im_shift_field_3, mask, C%mapped_im_field_3, lon_gcm, lat_gcm)
   ! Note that the im_field_* are here converted back to the IM units.
  END IF

  ! Finally the mapped fields gcm_field_* and the coordinates lon_gcm and lat_gcm are written to an output file:
  IF(C%choice_1D_lonlat_in_netcdf .EQV. .TRUE.) THEN
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: -
   CALL write_data_for_gcm(C%gcm_created_filename, lon_gcm, lat_gcm, C%mapped_gcm_field_1, gcm_field_1, &
                                                                     C%mapped_gcm_field_2, gcm_field_2, &
                                                                     C%mapped_gcm_field_3, gcm_field_3)
  ELSE
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: -
   CALL write_data_for_gcm_with_2D_lonlat(C%gcm_created_filename, lon_gcm, lat_gcm, C%mapped_gcm_field_1, gcm_field_1, &
                                                                                    C%mapped_gcm_field_2, gcm_field_2, &
                                                                                    C%mapped_gcm_field_3, gcm_field_3)
  END IF

  ! Finishing message:
  WRITE(UNIT=*, FMT='(3A)') ' Finished! The file  ', TRIM(C%gcm_created_filename), '  is created'

  CALL finalize_constants()
END PROGRAM oblimap_im_to_gcm_program



  SUBROUTINE initialize_extended_im_field(im_field, im_field_extended)
    ! To ensure that during the interpolation enough projected IM points surround the GCM grid 
    ! points, the IM grid is extended.  
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%NX,C%NY),                   INTENT(IN)  :: im_field

    ! Output variable:
    REAL(dp), DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED), INTENT(OUT) :: im_field_extended
        
    ! Local variables:
    INTEGER                                                       :: m, n

    ! Extending the IM domain with a margin (outside the C%NX-C%NY domain), and giving 
    ! the field at the margin the values of the nearest domain boundary point. 
      
    ! Filling the usual domain:
    im_field_extended(C%NX_ext+1:C%NX_ext+C%NX,C%NY_ext+1:C%NY_ext+C%NY) = im_field
    ! Filling the left margin:
    DO m = 1, C%NX_ext
     im_field_extended(m,C%NY_ext+1:C%NY_ext+C%NY) = im_field(1,:)
    END DO
    ! Filling the right margin:
    DO m = C%NX + C%NX_ext + 1, C%NX_EXTENDED
     im_field_extended(m,C%NY_ext+1:C%NY_ext+C%NY) = im_field(C%NX,:)
    END DO
    ! Filling the bottom margin:
    DO n = 1, C%NY_ext
     im_field_extended(:,n) = im_field_extended(:,C%NY_ext + 1)
    END DO
    ! Filling the upper margin:
    DO n = C%NY + C%NY_ext + 1, C%NY_EXTENDED
     im_field_extended(:,n) = im_field_extended(:,C%NY_ext + C%NY)
    END DO
  END SUBROUTINE initialize_extended_im_field
