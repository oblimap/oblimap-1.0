! File name: oblimap_gcm_to_im_program.f90          
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

PROGRAM oblimap_gcm_to_im_program
  ! This program reads several field variables from a Global Circulation Model (GCM) netcdf file,
  ! which are projected by an oblique stereographic projection on a rectangular Ice Model (IM)
  ! grid. The GCM grid point coordinates are projected on the IM grid, where in general they fall
  ! irregular and between the IM grid points. Therefore the field values defined on the GCM grid  
  ! points have to be interpolated at each IM grid point, with help of the nearby projected GCM
  ! grid points. Two interpolation methods are available in this program both based on a distance 
  ! weigthing Shepard technique. One method, the 'quadrant method', searches within each quadrant 
  ! around each IM grid point the nearest projected GCM point and interpolates with help of this
  ! four points to obtain the field value at such a IM grid point. Another method, the 'radius
  ! method', searches all projected GCM points within a certain radius and interpolates with help 
  ! of these points to obtain the field value at such a IM grid point. Finally the projected and
  ! interpolated IM fields are written to an IM netcdf file.
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
  USE oblimap_configuration_module, ONLY: dp, C, read_config_file, initialize_constants, finalize_constants
  USE oblimap_mapping_module, ONLY: fast_mapping
  USE oblimap_read_and_write_module, ONLY: read_data_gcm, read_data_gcm_with_2D_lonlat, write_data_for_im
  USE oblimap_scan_contributions_module, ONLY: scan_with_quadrant_method_gcm_to_im, scan_with_radius_method_gcm_to_im
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lon_gcm
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: lat_gcm
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_2
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gcm_field_3
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_2
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: im_field_3
  
  ! Read the config-variables which are declared in the grid module from the configuration file:
  CALL read_config_file()
  ! Initialization of the struckt  C% :
  CALL initialize_constants()

  ! Allocate all arrays for which the dimension sizes are given in the CONFIG file:
  ALLOCATE(lon_gcm            (C%NLON,C%NLAT))
  ALLOCATE(lat_gcm            (C%NLON,C%NLAT))
  ALLOCATE(gcm_field_1        (C%NLON,C%NLAT))
  ALLOCATE(gcm_field_2        (C%NLON,C%NLAT))
  ALLOCATE(gcm_field_3        (C%NLON,C%NLAT))
  ALLOCATE(im_field_1         (C%NX  ,C%NY  ))
  ALLOCATE(im_field_2         (C%NX  ,C%NY  ))
  ALLOCATE(im_field_3         (C%NX  ,C%NY  ))
  
  ! Inform message about the intersection angle alpha:
  IF(C%choice_projection_method == 'oblique_stereographic_projection' .OR. &
     C%choice_projection_method == 'oblique_stereographic_projection_snyder' .OR. &
     C%choice_projection_method == 'oblique_stereographic_projection_ellipsoid_snyder') THEN
   IF(REAL(C%NX * C%NY * C%dx * C%dy) > 2._dp * C%pi * C%R_earth**2) THEN
    WRITE(UNIT=*, FMT='((2A, F12.3))') &
     '\n No optimal alpha_stereographic_config could be calcultated for this grid configuration. The best estimate might \n', &
     ' be alpha_stereographic_config = 90 degree, or take a smaller grid configuration. This run uses alpha = ', C%rad2deg * C%alpha_stereographic, '\n'
   ELSE
    WRITE(UNIT=*, FMT='(2(A, F12.3))') &
     'An optimal alpha_stereographic_config = ', C%rad2deg * ASIN(SQRT((C%NX * C%dx * C%NY * C%dy) / (2._dp * C%pi)) / C%R_earth), &
     ' degree. This run uses alpha = ', C%rad2deg * C%alpha_stereographic
   END IF
  END IF

  ! Reading the GCM field variables 'C%mapped_gcm_field_*' and the longitude and latitude coordinates
  ! of the GCM grid from a GCM-netcdf file.
  IF(C%choice_1D_lonlat_in_netcdf .EQV. .TRUE.) THEN
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: lon_gcm, lat_gcm, gcm_field_1, gcm_field_2, gcm_field_3
   CALL read_data_gcm(C%gcm_input_filename, C%mapped_gcm_field_1, lon_gcm, lat_gcm, gcm_field_1, &
                                            C%mapped_gcm_field_2,                   gcm_field_2, &
                                            C%mapped_gcm_field_3,                   gcm_field_3)
                                            
  ELSE
   ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
   ! Output: lon_gcm, lat_gcm, gcm_field_1, gcm_field_2, gcm_field_3
   CALL read_data_gcm_with_2D_lonlat(C%gcm_input_filename, C%mapped_gcm_field_1, lon_gcm, lat_gcm, gcm_field_1, &
                                                           C%mapped_gcm_field_2,                   gcm_field_2, &
                                                           C%mapped_gcm_field_3,                   gcm_field_3)

   ! Adhoc RACMO coordinate shift correction:
   WHERE(lon_gcm < 0._dp) lon_gcm = lon_gcm + 360._dp
  END IF
  
  ! Before the fast mapping routine can be used, the C%input_fast_map_GCM_to_IM file has to be created. This file
  ! contains for each target grid point the coordinates and the relative distances of the contributing points which 
  ! are used for the field interpolation. This contributing points of one such a target point are the nearest
  ! projected points selected by the radius method or by the quadrant method. This 'scan' part in fact contains the
  ! most technical and CPU consuming part of OBLIMAP covering the projection and the selection of points for the 
  ! interpolation (and keeping the relative distances of each projected point relative to the target point).
  IF(C%choice_fast_map_GGM_to_IM .EQV. .FALSE.) THEN
   IF(C%choice_radius_GCM_to_IM .EQV. .TRUE.) THEN
    ! Output: the C%input_fast_map_GCM_to_IM is created
    CALL scan_with_radius_method_gcm_to_im(lon_gcm, lat_gcm)
   ELSE
    ! Output: the C%input_fast_map_GCM_to_IM is created
    CALL scan_with_quadrant_method_gcm_to_im(lon_gcm, lat_gcm)
   END IF                                 
  END IF
  
  ! The GCM fields are mapped (= projected + interpolated) on the IM grid. For each target grid point the coordinates and 
  ! the relative distances of the nearest projected points are stored in the C%input_fast_map_GCM_to_IM file, these are 
  ! used to map the fields with the fast mapping:
  ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
  ! Output: im_field_1, im_field_2, im_field_3
  CALL fast_mapping(C%input_fast_map_GCM_to_IM, C%NLON, C%NLAT, C%NX, C%NY, gcm_field_1, im_field_1, &
                                                                            gcm_field_2, im_field_2, &
                                                                            gcm_field_3, im_field_3)

  ! Rescaling the fields by multiplication with a gcm_to_im_factor and by a gcm_to_im_shift before writing (in case the units differ):
  im_field_1 = C%gcm_to_im_factor_field_1 * im_field_1 + C%gcm_to_im_shift_field_1
  im_field_2 = C%gcm_to_im_factor_field_2 * im_field_2 + C%gcm_to_im_shift_field_2
  im_field_3 = C%gcm_to_im_factor_field_3 * im_field_3 + C%gcm_to_im_shift_field_3

  ! Finally the mapped fields im_field_* are written to an output file:
  ! Note that currently up to *_field_9 OPTIONAL arguments are possible.
  ! Output: -
  CALL write_data_for_im(C%im_created_filename, C%mapped_im_field_1, im_field_1, &
                                                C%mapped_im_field_2, im_field_2, &
                                                C%mapped_im_field_3, im_field_3)

  ! Finishing message:
  WRITE(UNIT=*, FMT='(3A)') ' Finished! The file  ', TRIM(C%im_created_filename), '  is created'

  CALL finalize_constants()
END PROGRAM oblimap_gcm_to_im_program

