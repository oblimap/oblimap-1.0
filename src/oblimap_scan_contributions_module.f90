! File name: oblimap_scan_contributions_module.f90          
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

MODULE oblimap_scan_contributions_module
  USE oblimap_configuration_module, ONLY: dp
  IMPLICIT NONE

  TYPE triplet
    INTEGER  :: row_index     ! row index of nearest point within one quadrant
    INTEGER  :: column_index  ! column index of nearest point within one quadrant
    REAL(dp) :: distance      ! distance of this nearest point relative to the IM point (m,n)
  END TYPE triplet
  


CONTAINS

  ! -----------------------------------------------------------------------------
  ! ROUTINES WHICH SCAN THE CONTRIBUTING POINTS FOR INTERPOLATION OF GCM TO IM
  ! -----------------------------------------------------------------------------

  SUBROUTINE scan_with_quadrant_method_gcm_to_im(lon_gcm, lat_gcm)
    ! This routine selects the contributing points for each target grid point, by searching with the quadrant method. First
    ! the coordinates of the GCM grid points are projected with the oblique stereographic projection to the IM coordinates. 
    ! Thereafter with these projected coordinates the distances of the projected points relative to each target grid point
    ! are calculated and used to select the nearest contributing grid points. The GCM-grid indices of the contributing points
    ! and the relative distance to 'their' target grid point are stored by writing them to the C%input_fast_map_GCM_to_IM 
    ! file. With the indices and the distances of the contributing points the GCM fields can be mapped fast and simultaneously
    ! on to the IM grid. 
    USE oblimap_configuration_module, ONLY : dp, C 
    IMPLICIT NONE

    ! Input variables:
    REAL(dp),      DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm                          ! longitude coordinates (degrees) of GCM grid
    REAL(dp),      DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm                          ! latitude coordinates  (degrees) of GCM grid

    ! Local variables:
    REAL(dp),      DIMENSION(C%NLON,C%NLAT)              :: x_coordinates_of_gcm_grid_points ! The x-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(C%NLON,C%NLAT)              :: y_coordinates_of_gcm_grid_points ! The y-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(C%NX,C%NY)                  :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(C%NX,C%NY)                  :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'
    INTEGER                                              :: i, j
    INTEGER                                              :: m, n
    INTEGER                                              :: quadrant                         ! The quadrant I, II, III or IV relative to an IM grid point
    TYPE(triplet)                                        :: projected_gcm                    ! Projected GCM point on S' 
    TYPE(triplet), DIMENSION(4)                          :: contribution                     ! Nearest projected GCM point in quadrant (DIM=I,II,III or IV) in S', relative to the IM grid point
    
    ! Projection of the GCM coordinates to the IM coordinates with the oblique stereographic projection:
    ! Output: x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points
    CALL gcm_projected_xy_coordinates(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)
    
    ! Opening the 'fast input'  file:
    OPEN(UNIT=107, FILE=TRIM(C%input_fast_map_GCM_to_IM))
    
    ! Writing the header of the C%input_fast_map_GCM_to_IM file:
    WRITE(UNIT=107, FMT='(A)') '# Do not remove this header. The quadrant method was used. The format is:'
    WRITE(UNIT=107, FMT='(A)') '#  m  n  N  N(i  j  distance)'
    WRITE(UNIT=107, FMT='(A)') '# with i = the longitudinal GCM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with j = the latitudinal GCM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with N = amount of weighted poits'
    WRITE(UNIT=107, FMT='(A)') '# with m = the x-axis IM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with n = the y-axis IM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# and distance is the distance between the GCM and IM points'
    WRITE(UNIT=107, FMT='(A, 3(A, F8.4), 2(A, I6))')            '# ', ' alpha_stereographic = ', C%rad2deg * C%alpha_stereographic, ', lambda_M = ', C%rad2deg * C%lambda_M, ', phi_M = ', C%rad2deg * C%phi_M, ', NLON = ', C%NLON, ', NLAT = ', C%NLAT
    WRITE(UNIT=107, FMT='(A, 2(A, I6), 2(A, F12.2), 2(A, I6))') '# ', ' NX = ', C%NX, ', NY = ', C%NY, ', dx = ', C%dx, ', dy = ', C%dy, ', NX_ext = ', C%NX_ext, ', NY_ext = ', C%NY_ext
    WRITE(UNIT=107, FMT='(A)') '# The first number below is the amount of mapped GCM points:'

    ! Number of lines in this the 'fast input'  file, equal to the amount mapped GCM points:
    WRITE(UNIT=107, FMT='(I10)') C%NX * C%NY 

    ! Determine the x, y coordinates of the IM grid points in plane S'
    ! Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    CALL initialize_map_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points) 

    ! For each IM grid point the four closest projected GCM points are determined:
    WRITE(UNIT=*,FMT='(A,I3,A)') '  m =   1 (last m  = ',C%NX,')'
    DO m = 1, C%NX
      IF(MOD(m*1.0,10.0) == 0) WRITE(UNIT=*,FMT='(A,I3)') '  m = ', m
    DO n = 1, C%NY
    
      ! Initialize the four nearest quadrant distances to a large value (the width of the IM grid domain):
      contribution(:)%distance = C%dx * C%NX
      
      DO i = 1, C%NLON
      DO j = 1, C%NLAT
        ! Determine the quadrant in which the projected point lies relative to the considered grid point: 
        ! Output: quadrant
        CALL find_quadrant_around_IM_grid_point(x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j), &
                                                x_coordinates_of_im_grid_points(m,n),  y_coordinates_of_im_grid_points(m,n), quadrant)
        
        ! Determine in the flat plane S' the distance between the projected GCM coordinates relative to the considered IM grid point: 
        projected_gcm%row_index    = i
        projected_gcm%column_index = j
        projected_gcm%distance     = SQRT((x_coordinates_of_gcm_grid_points(i,j) - x_coordinates_of_im_grid_points(m,n))**2 + &
                                          (y_coordinates_of_gcm_grid_points(i,j) - y_coordinates_of_im_grid_points(m,n))**2)
        ! In case the projected point coincides with the grid point we put it at the very close distance of 1 centimeter, preventing devision by zero:
        IF(projected_gcm%distance == 0._dp) projected_gcm%distance = 0.01_dp
        
        ! Select the in S' projected GCM point with the shortest distance to the considered IM grid point in this quadrant, 
        ! and keep this distance and the GCM-grid indices of this GCM point in S:
        IF(projected_gcm%distance < contribution(quadrant)%distance) contribution(quadrant) = projected_gcm
      END DO
      END DO

      WRITE(UNIT=107,FMT='(3I6, 4(2I6,E23.15))') m, n, 4, &
       contribution(1)%row_index, contribution(1)%column_index, contribution(1)%distance, &
       contribution(2)%row_index, contribution(2)%column_index, contribution(2)%distance, &
       contribution(3)%row_index, contribution(3)%column_index, contribution(3)%distance, &
       contribution(4)%row_index, contribution(4)%column_index, contribution(4)%distance
    END DO
    END DO

    ! Closing the the 'fast input'  file:
    CLOSE(UNIT=107)
  END SUBROUTINE scan_with_quadrant_method_gcm_to_im



  SUBROUTINE gcm_projected_xy_coordinates(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)
    ! This routine projects the GCM coordinates on the requested plane S' which coincides with the IM grid,
    ! with an oblique stereographic projection. 
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_sg_projection_ellipsoid_snyder, &
      oblique_laea_projection_snyder, oblique_laea_projection_ellipsoid_snyder
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm                          ! longitude coordinates (degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm                          ! latitude coordinates  (degrees) of GCM grid

    ! Output variable:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: x_coordinates_of_gcm_grid_points ! The x-coordinates of the GCM points projected om S'
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: y_coordinates_of_gcm_grid_points ! The y-coordinates of the GCM points projected om S'

    ! Local variables:
    INTEGER                                         :: i, j
    
    ! Determine the x,y coordinates of each GCM longitude-latitude coordinate after 
    ! The oblique stereographic projection on the projection plane S'
    DO i = 1, C%NLON
    DO j = 1, C%NLAT
      IF(C%choice_projection_method == 'oblique_stereographic_projection' .OR. &
         C%choice_projection_method == 'oblique_stereographic_projection_snyder') THEN
       ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
       CALL oblique_sg_projection(lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      ELSE IF(C%choice_projection_method == 'oblique_stereographic_projection_ellipsoid_snyder') THEN
       ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
       CALL oblique_sg_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_snyder') THEN
       ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
       CALL oblique_laea_projection_snyder(lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_ellipsoid_snyder') THEN
       ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
       CALL oblique_laea_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      ELSE
       STOP ' Wrong choice in CONFIG file for C%choice_projection_method [gcm_projected_xy_coordinates], --STOPPED'
      END IF
    END DO
    END DO
  END SUBROUTINE gcm_projected_xy_coordinates



  SUBROUTINE initialize_map_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)
    ! Initialize the x, y coordinate values of the IM grid points in plane S'. In total there are
    ! C%NX * C%NY grid points and the grid cell distances are C%dx and C%dy (in meters). The central 
    ! grid point is (0,0) and coincide with the longitude-latitude coordinates (lamda_M, phi_M).
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Output variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: x_coordinates_of_im_grid_points    
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: y_coordinates_of_im_grid_points    

    ! Local variables:
    INTEGER                                     :: m, n

    DO m = 1, C%NX
    DO n = 1, C%NY
      x_coordinates_of_im_grid_points(m,n) = C%dx * (m - ((C%NX+1) / 2))
      y_coordinates_of_im_grid_points(m,n) = C%dy * (n - ((C%NY+1) / 2))
    END DO
    END DO
  END SUBROUTINE initialize_map_gcm_to_im



  SUBROUTINE find_quadrant_around_IM_grid_point(x_value_projected_point,  y_value_projected_point, &
                                                x_value_IM_grid_point, y_value_IM_grid_point, quadrant)              
    ! Determing the quadrant in which a 'projected point' is situated relative to an IM grid point.                                           
    !   quadrants:
    !   II  |   I
    !       |     
    !  -----|-----
    !       |     
    !  III  |  IV 
    !
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE                                                                                      
                              
    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_value_projected_point
    REAL(dp), INTENT(IN)  :: y_value_projected_point
    REAL(dp), INTENT(IN)  :: x_value_IM_grid_point
    REAL(dp), INTENT(IN)  :: y_value_IM_grid_point

    ! Output variable:
    INTEGER,  INTENT(OUT) :: quadrant 

    IF(      x_value_projected_point >  x_value_IM_grid_point .AND. y_value_projected_point >= y_value_IM_grid_point) THEN
     ! Check for quadrant I
     quadrant = 1
    ELSE IF( x_value_projected_point <= x_value_IM_grid_point .AND. y_value_projected_point >  y_value_IM_grid_point) THEN
     ! Check for quadrant II
     quadrant = 2
    ELSE IF( x_value_projected_point <  x_value_IM_grid_point .AND. y_value_projected_point <= y_value_IM_grid_point) THEN
     ! Check for quadrant III
     quadrant = 3
    ELSE IF( x_value_projected_point >= x_value_IM_grid_point .AND. y_value_projected_point <= y_value_IM_grid_point) THEN
     ! Check for quadrant IV
     quadrant = 4
    END IF
  END SUBROUTINE find_quadrant_around_IM_grid_point
  
  
  
  SUBROUTINE scan_with_radius_method_gcm_to_im(lon_gcm, lat_gcm)
    ! This routine selects the contributing points for each target grid point, by searching with the radius method. First
    ! the coordinates of the GCM grid points are projected with the oblique stereographic projection to the IM coordinates. 
    ! Thereafter with these projected coordinates the distances of the projected points relative to each target grid point
    ! are calculated and used to select the nearest contributing grid points. The GCM-grid indices of the contributing points
    ! and the relative distance to 'their' target grid point are stored by writing them to the C%input_fast_map_GCM_to_IM 
    ! file. With the indices and the distances of the contributing points the GCM fields can be mapped fast and simultaneously
    ! on to the IM grid. 
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE             

    ! Input variables:
    REAL(dp),      DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm                          ! longitude coordinates (degrees) of GCM grid
    REAL(dp),      DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm                          ! latitude coordinates  (degrees) of GCM grid

    ! Local variables:
    INTEGER                                              :: max_size
    REAL(dp),      DIMENSION(C%NLON,C%NLAT)              :: x_coordinates_of_gcm_grid_points ! The x-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(C%NLON,C%NLAT)              :: y_coordinates_of_gcm_grid_points ! The y-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(C%NX,C%NY)                  :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(C%NX,C%NY)                  :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'
    INTEGER                                              :: i, j
    INTEGER                                              :: m, n
    INTEGER                                              :: count_contributions
    INTEGER                                              :: maximum_contributions = 0
    INTEGER                                              :: loop
    TYPE(triplet)                                        :: projected_gcm                    ! Projected GCM point on S' 
    TYPE(triplet), DIMENSION(:), ALLOCATABLE             :: contribution                     ! Nearest projected GCM point in quadrant (DIM=I,II,III or IV) in S', relative to the IM grid point
    
    ! The devision by 1000 is to prevent the failure of CEILING with latge numbers:
    max_size = CEILING(MAX(4._dp * C%pi * (C%R_search_GCM_to_IM / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)), &
                                   C%pi * (C%R_search_GCM_to_IM / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)))  * C%oblimap_allocate_factor, dp)
    ALLOCATE(contribution(max_size))
    WRITE(*,*) 'max_size = ', max_size
    
    ! Projection of the GCM coordinates to the IM coordinates with the oblique stereographic projection:
    ! Output: x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points
    CALL gcm_projected_xy_coordinates(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)
    
    ! Opening the 'fast input'  file:
    OPEN(UNIT=107, FILE=TRIM(C%input_fast_map_GCM_to_IM))
    
    ! Writing the header of the C%input_fast_map_GCM_to_IM file:
    WRITE(UNIT=107, FMT='(A)') '# Do not remove this header. The radius method was used. The format is:'
    WRITE(UNIT=107, FMT='(A)') '#  m  n  N  N(i  j  distance)'
    WRITE(UNIT=107, FMT='(A)') '# with i = the longitudinal GCM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with j = the latitudinal GCM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with N = amount of weighted poits'
    WRITE(UNIT=107, FMT='(A)') '# with m = the x-axis IM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# with n = the y-axis IM grid counter'
    WRITE(UNIT=107, FMT='(A)') '# and distance is the distance between the GCM and IM points'
    WRITE(UNIT=107, FMT='(A, 3(A, F8.4), 2(A, I6))')            '# ', ' alpha_stereographic = ', C%rad2deg * C%alpha_stereographic, ', lambda_M = ', C%rad2deg * C%lambda_M, ', phi_M = ', C%rad2deg * C%phi_M, ', NLON = ', C%NLON, ', NLAT = ', C%NLAT
    WRITE(UNIT=107, FMT='(A, 2(A, I6), 2(A, F12.2), 2(A, I6))') '# ', ' NX = ', C%NX, ', NY = ', C%NY, ', dx = ', C%dx, ', dy = ', C%dy, ', NX_ext = ', C%NX_ext, ', NY_ext = ', C%NY_ext
    WRITE(UNIT=107, FMT='(A)') '# The first number below is the amount of mapped GCM points:'

    ! Number of lines in this the 'fast input'  file, equal to the amount mapped GCM points:
    WRITE(UNIT=107, FMT='(I10)') C%NX * C%NY 

    ! Determine the x, y coordinates of the IM grid points in plane S'
    ! Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    CALL initialize_map_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points) 

    ! For each IM grid point the four closest projected GCM points are determined:
    WRITE(UNIT=*,FMT='(A,I3,A)') '  m =   1 (last m  = ',C%NX,')'
    DO m = 1, C%NX
      IF(MOD(m*1.0,10.0) == 0) WRITE(UNIT=*,FMT='(A,I3)') '  m = ', m
    DO n = 1, C%NY

      ! See equation (2.19) in Reerink et al. (2010):  
      count_contributions = 0
      DO i = 1, C%NLON
      DO j = 1, C%NLAT
         projected_gcm%row_index    = i
         projected_gcm%column_index = j
         projected_gcm%distance     = SQRT((x_coordinates_of_gcm_grid_points(i,j) - x_coordinates_of_im_grid_points(m,n))**2 + &
                                           (y_coordinates_of_gcm_grid_points(i,j) - y_coordinates_of_im_grid_points(m,n))**2)

         IF(projected_gcm%distance <= C%R_search_GCM_to_IM .AND. projected_gcm%distance > 0._dp) THEN
          count_contributions = count_contributions + 1 
          maximum_contributions = MAX(maximum_contributions, count_contributions)
          IF(count_contributions > max_size) WRITE(UNIT=*,FMT='(2(A, I10))') 'WARNING: the array contribution is not allocated properly [scan_with_radius_method_gcm_to_im], number of contributions =', count_contributions, ', max size = ', max_size
          contribution(count_contributions) = projected_gcm
         END IF
      END DO
      END DO
      WRITE(UNIT=107, FMT='(3I6)', ADVANCE='NO') m, n, count_contributions
      DO loop = 1, count_contributions
         WRITE(UNIT=107, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop)%row_index, contribution(loop)%column_index, contribution(loop)%distance
      END DO
      WRITE(UNIT=107, FMT='(A)') ''

      IF(count_contributions == 0) &
       WRITE(UNIT=*,FMT='(A, I10)') 'WARNING: There are 0 points within C%R_search_GCM_to_IM = ', C%R_search_GCM_to_IM
    END DO
    END DO    
    WRITE(UNIT=*, FMT='(A, I5)') ' Maximum number of weighted points per IM grid point = ', maximum_contributions
    IF(maximum_contributions > max_size) THEN
     WRITE(UNIT=*, FMT='(A, F4.1)') '\n OBLIMAP has been STOPPED: The oblimap_allocate_factor_config should be increased to ', 1.1_dp * maximum_contributions / REAL(max_size), 'in the CONFIG file.\n'
     STOP ''
    END IF

    ! Closing the the 'fast input'  file:
    CLOSE(UNIT=107)
  END SUBROUTINE scan_with_radius_method_gcm_to_im



  ! -----------------------------------------------------------------------------
  ! ROUTINES WHICH SCAN THE CONTRIBUTING POINTS FOR INTERPOLATION OF IM TO GCM
  ! -----------------------------------------------------------------------------


  SUBROUTINE scan_with_radius_method_im_to_gcm(lon_gcm, lat_gcm)
    ! This routine selects the contributing points for each target grid point, by searching with the radius method. First
    ! the coordinates of the IM grid points are projected with the inverse oblique stereographic projection to the GCM 
    ! coordinates. Thereafter with these projected coordinates the distances of the projected points relative to each target
    ! grid point are calculated and used to select the nearest contributing grid points. The IM-grid indices of the
    ! contributing points and the relative distance to 'their' target grid point are stored by writing them to the  
    ! C%input_fast_map_IM_to_GCM file. With the indices and the distances of the contributing points the IM fields can be
    ! mapped fast and simultaneously on to the GCM grid. 
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp),      DIMENSION(C%NLON,C%NLAT),               INTENT(IN)  :: lon_gcm
    REAL(dp),      DIMENSION(C%NLON,C%NLAT),               INTENT(IN)  :: lat_gcm
           
    ! Local variables:
    INTEGER                                                            :: max_size
    REAL(dp),      DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED)              :: lon_coordinates_of_im_grid_points
    REAL(dp),      DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED)              :: lat_coordinates_of_im_grid_points
    INTEGER,       DIMENSION(C%NLON,C%NLAT)                            :: mask 
    INTEGER                                                            :: i, j
    INTEGER                                                            :: m, n
    INTEGER                                                            :: count_contributions
    INTEGER                                                            :: maximum_contributions = 0
    INTEGER                                                            :: loop
    TYPE(triplet)                                                      :: projected_im   ! Projected IM point on S 
    TYPE(triplet), DIMENSION(:), ALLOCATABLE                           :: contribution   ! Projected IM points on S within C%R_search_IM_to_GCM, relative to the GCM grid point
    
    max_size = CEILING(MAX(4._dp * C%pi * (C%R_search_IM_to_GCM / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)), &
                                   C%pi * (C%R_search_IM_to_GCM / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)))  * C%oblimap_allocate_factor, dp)
    ALLOCATE(contribution(max_size))

    ! Projection of the IM coordinates to the GCM coordinates with the inverse oblique stereographic projection:
    ! Output: lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points
    CALL im_projected_lonlat_coordinates(lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    
    ! A mask is created for the GCM points which are involved (mask ==1) in the mapping:
    ! Output: mask
    CALL make_mask(lon_gcm, lat_gcm, mask)
    
    ! Opening the 'fast input' file:
    OPEN(UNIT=108, FILE=TRIM(C%input_fast_map_IM_to_GCM))
    
    ! Writing the header of the C%input_fast_map_GCM_to_IM file:
    WRITE(UNIT=108, FMT='(A)') '# Do not remove this header. The radius method was used. The format is:'
    WRITE(UNIT=108, FMT='(A)') '#  i  j  N  N(m  n  distance)'
    WRITE(UNIT=108, FMT='(A)') '# with i = the longitudinal GCM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with j = the latitudinal GCM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with N = amount of weighted poits'
    WRITE(UNIT=108, FMT='(A)') '# with m = the x-axis IM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with n = the y-axis IM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# and distance is the distance between the GCM and IM points'
    WRITE(UNIT=108, FMT='(A, 3(A, F8.4), 2(A, I6))')            '# ', ' alpha_stereographic = ', C%rad2deg * C%alpha_stereographic, ', lambda_M = ', C%rad2deg * C%lambda_M, ', phi_M = ', C%rad2deg * C%phi_M, ', NLON = ', C%NLON, ', NLAT = ', C%NLAT
    WRITE(UNIT=108, FMT='(A, 2(A, I6), 2(A, F12.2), 2(A, I6))') '# ', ' NX = ', C%NX, ', NY = ', C%NY, ', dx = ', C%dx, ', dy = ', C%dy, ', NX_ext = ', C%NX_ext, ', NY_ext = ', C%NY_ext
    WRITE(UNIT=108, FMT='(A)') '# The first number below is the amount of mapped GCM points:'

    ! Number of lines in this 'fast input'  file, equal to the amount mapped IM points:
    WRITE(UNIT=108, FMT='(I10)') COUNT(mask==1)
    
    WRITE(UNIT=*,FMT='(A,I3,A)') '  i =   1 (last i  = ',C%NLON,')'
    DO i = 1, C%NLON
      IF(MOD(i*1.0,10.0) == 0)  WRITE(UNIT=*,FMT='(A,I3)') '  i = ', i
    DO j = 1, C%NLAT
      
      ! See equation (2.19) in Reerink et al. (2010):  
      count_contributions = 0 
      IF(mask(i,j) == 1) THEN
       DO m = 1, C%NX_EXTENDED
       DO n = 1, C%NY_EXTENDED
         projected_im%row_index    = m
         projected_im%column_index = n
         projected_im%distance     = distance_over_curved_surface(lon_gcm(i,j), lat_gcm(i,j), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
         IF(projected_im%distance <= C%R_search_IM_to_GCM .AND. projected_im%distance > 0._dp) THEN
          count_contributions = count_contributions + 1 
          maximum_contributions = MAX(maximum_contributions, count_contributions)
          IF(count_contributions > max_size) WRITE(UNIT=*,FMT='(2(A, I10))') 'WARNING: the array contribution is not allocated properly [scan_with_radius_method_im_to_gcm], number of contributions =', count_contributions, ', max size = ', max_size
          contribution(count_contributions) = projected_im
         END IF
       END DO
       END DO
       WRITE(UNIT=108, FMT='(3I6)', ADVANCE='NO') i, j, count_contributions
       DO loop = 1, count_contributions
          WRITE(UNIT=108, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop)%row_index, contribution(loop)%column_index, contribution(loop)%distance
       END DO
       WRITE(UNIT=108, FMT='(A)') ''

       IF(count_contributions == 0) WRITE(UNIT=*,FMT='(A)') 'WARNING: There are 0 points within C%R_search_IM_to_GCM'
      END IF
    END DO
    END DO  
    WRITE(UNIT=*, FMT='(A, I5)') ' Maximum number of weighted points per GCM grid point = ', maximum_contributions
    IF(maximum_contributions > max_size) THEN
     WRITE(UNIT=*, FMT='(A, F4.1)') '\n OBLIMAP has been STOPPED: The oblimap_allocate_factor_config should be increased to ', 1.1_dp * maximum_contributions / REAL(max_size), 'in the CONFIG file.\n'
     STOP ''
    END IF
    
    ! Closing the 'fast input' file:
    CLOSE(UNIT=108)
  END SUBROUTINE scan_with_radius_method_im_to_gcm  



  SUBROUTINE im_projected_lonlat_coordinates(lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    ! This routine projects the IM coordinates on the requested plane S'which coincides with the GCM grid,
    ! with an inverse oblique stereographic projection. 
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: inverse_oblique_sg_projection, inverse_oblique_sg_projection_snyder, &
      inverse_oblique_sg_projection_ellipsoid_snyder, &
      inverse_oblique_laea_projection_snyder, inverse_oblique_laea_projection_ellipsoid_snyder
    IMPLICIT NONE

    ! Output variable:
    REAL(dp), DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED), INTENT(OUT) :: lon_coordinates_of_im_grid_points
    REAL(dp), DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED), INTENT(OUT) :: lat_coordinates_of_im_grid_points

    ! Local variables:
    INTEGER                                                       :: m, n
    REAL(dp), DIMENSION(C%NX_EXTENDED)                            :: x_coordinates_of_im_extended_grid_points
    REAL(dp), DIMENSION(              C%NY_EXTENDED)              :: y_coordinates_of_im_extended_grid_points

    ! Output: x_coordinates_of_im_extended_grid_points, y_coordinates_of_im_extended_grid_points
    CALL initialize_extend_im_coordinates(x_coordinates_of_im_extended_grid_points, y_coordinates_of_im_extended_grid_points)
     
    ! Project the (extended) IM (x,y) -coordinates to the longitude-latitude values:  
    DO m = 1, C%NX_EXTENDED
    DO n = 1, C%NY_EXTENDED
      IF(C%choice_projection_method == 'oblique_stereographic_projection') THEN
       ! Output: im_ext_grid_point_lon_coordinate(m,n), im_ext_grid_point_lat_coordinate(m,n) 
       CALL inverse_oblique_sg_projection(x_coordinates_of_im_extended_grid_points(m), y_coordinates_of_im_extended_grid_points(n), &
                                          lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      ELSE IF(C%choice_projection_method == 'oblique_stereographic_projection_snyder') THEN
       ! Output: im_ext_grid_point_lon_coordinate(m,n), im_ext_grid_point_lat_coordinate(m,n) 
       CALL inverse_oblique_sg_projection_snyder(x_coordinates_of_im_extended_grid_points(m), y_coordinates_of_im_extended_grid_points(n), &
                                                 lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      ELSE IF(C%choice_projection_method == 'oblique_stereographic_projection_ellipsoid_snyder') THEN
       ! Output: im_ext_grid_point_lon_coordinate(m,n), im_ext_grid_point_lat_coordinate(m,n) 
       CALL inverse_oblique_sg_projection_ellipsoid_snyder(x_coordinates_of_im_extended_grid_points(m), y_coordinates_of_im_extended_grid_points(n), &
                                                            lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_snyder') THEN
       ! Output: im_ext_grid_point_lon_coordinate(m,n), im_ext_grid_point_lat_coordinate(m,n) 
       CALL inverse_oblique_laea_projection_snyder(x_coordinates_of_im_extended_grid_points(m), y_coordinates_of_im_extended_grid_points(n), &
                                                   lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_ellipsoid_snyder') THEN
       ! Output: im_ext_grid_point_lon_coordinate(m,n), im_ext_grid_point_lat_coordinate(m,n) 
       CALL inverse_oblique_laea_projection_ellipsoid_snyder(x_coordinates_of_im_extended_grid_points(m), y_coordinates_of_im_extended_grid_points(n), &
                                                              lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      ELSE
       STOP ' Wrong choice in CONFIG file for C%choice_projection_method [im_projected_lonlat_coordinates], --STOPPED'
      END IF
    END DO
    END DO
  END SUBROUTINE im_projected_lonlat_coordinates



  SUBROUTINE initialize_extend_im_coordinates(x_coordinates_of_im_extended_grid_points, y_coordinates_of_im_extended_grid_points)
    ! To ensure that during the interpolation enough projected IM points surround the GCM grid 
    ! points, the IM grid is extended.  
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Output variable:
    REAL(dp), DIMENSION(C%NX_EXTENDED), INTENT(OUT) :: x_coordinates_of_im_extended_grid_points
    REAL(dp), DIMENSION(C%NY_EXTENDED), INTENT(OUT) :: y_coordinates_of_im_extended_grid_points
        
    ! Local variables:
    INTEGER                                         :: m, n
    REAL(dp)                                        :: extended_boundary_x
    REAL(dp)                                        :: extended_boundary_y

    extended_boundary_x = - ((C%NX - 1) / 2 + C%NX_ext) * C%dx
    extended_boundary_y = - ((C%NY - 1) / 2 + C%NY_ext) * C%dy

    ! Determine the x coordinate values of the extended grid (in meters).  
    DO m = 1, C%NX_EXTENDED
      x_coordinates_of_im_extended_grid_points(m) = extended_boundary_x + C%dx * (m - 1)
    END DO
    ! Determine the y coordinate values of the extended grid (in meters).  
    DO n = 1, C%NY_EXTENDED
      y_coordinates_of_im_extended_grid_points(n) = extended_boundary_y + C%dy * (n - 1)
    END DO
  END SUBROUTINE initialize_extend_im_coordinates
   
   

  SUBROUTINE make_mask(lon_gcm, lat_gcm, mask)  
    ! This routine determines the GCM longitude, latitude grid points that participate
    ! in the mapping. It distinguishes with a mask between points that are projected within the
    ! IM grid domain and points which are projected outside that domain.
    !
    ! With the mapping from the IM to the GCM (the inverse oblique stereographic projection 
    ! followed by the interpolation) it is necessary to know which GCM points are affected 
    ! by the mapping, the rest of the GCM points remain their initial value
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_sg_projection_ellipsoid_snyder, &
      oblique_laea_projection_snyder, oblique_laea_projection_ellipsoid_snyder
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm

    ! Output variable
    INTEGER,  DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: mask

    ! Local variables:
    INTEGER                                         :: i, j
    REAL(dp)                                        :: x_gcm
    REAL(dp)                                        :: y_gcm
    REAL(dp)                                        :: boundary_x
    REAL(dp)                                        :: boundary_y
 
    ! Determine the x and y coordinates of the boundaries of the original IM grid (in meters).  
    boundary_x  = ((C%NX - 1) / 2) * C%dx
    boundary_y  = ((C%NY - 1) / 2) * C%dy
 
    mask = 0
    DO i = 1, C%NLON
    DO j = 1, C%NLAT
      IF(C%choice_projection_method == 'oblique_stereographic_projection' .OR. &
         C%choice_projection_method == 'oblique_stereographic_projection_snyder') THEN
       ! Output: x_gcm, y_gcm
       CALL oblique_sg_projection(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      ELSE IF(C%choice_projection_method == 'oblique_stereographic_projection_ellipsoid_snyder') THEN
       ! Output: x_gcm, y_gcm
       CALL oblique_sg_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_snyder') THEN
       ! Output: x_gcm, y_gcm
       CALL oblique_laea_projection_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      ELSE IF(C%choice_projection_method == 'oblique_lambert_equal-area_projection_ellipsoid_snyder') THEN
       ! Output: x_gcm, y_gcm
       CALL oblique_laea_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      ELSE
       STOP ' Wrong choice in CONFIG file for C%choice_projection_method [make_mask], --STOPPED'
      END IF
      IF(x_gcm < boundary_x .AND. x_gcm > - boundary_x .AND. y_gcm < boundary_y .AND. y_gcm > - boundary_y) &
       mask(i,j) = 1
    END DO
    END DO
  END SUBROUTINE make_mask



  REAL(dp) FUNCTION  distance_over_curved_surface(point_1_lon, point_1_lat, point_2_lon, point_2_lat) RESULT (distance)
    ! Calculation of the distance over the Earth surface of two points which are situated at this earth surface.
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) ::point_1_lon 
    REAL(dp), INTENT(IN) ::point_1_lat 
    REAL(dp), INTENT(IN) ::point_2_lon 
    REAL(dp), INTENT(IN) ::point_2_lat 

   ! See equation (2.18) in Reerink et al. (2010):
   distance = C%R_earth * ACOS(COS(C%deg2rad * point_1_lat) * COS(C%deg2rad * point_2_lat) * COS(C%deg2rad * (point_1_lon - point_2_lon)) + &
                               SIN(C%deg2rad * point_1_lat) * SIN(C%deg2rad * point_2_lat))
  END FUNCTION distance_over_curved_surface
   


  SUBROUTINE scan_with_quadrant_method_im_to_gcm(lon_gcm, lat_gcm)
    ! This routine selects the contributing points for each target grid point, by searching with the quadrant method. First
    ! the coordinates of the IM grid points are projected with the inverse oblique stereographic projection to the GCM 
    ! coordinates. Thereafter with these projected coordinates the distances of the projected points relative to each target
    ! grid point are calculated and used to select the nearest contributing grid points. The IM-grid indices of the
    ! contributing points and the relative distance to 'their' target grid point are stored by writing them to the  
    ! C%input_fast_map_IM_to_GCM file. With the indices and the distances of the contributing points the IM fields can be
    ! mapped fast and simultaneously on to the GCM grid. 
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE
    
    ! Input variables:
    REAL(dp),      DIMENSION(C%NLON,C%NLAT),               INTENT(IN)  :: lon_gcm
    REAL(dp),      DIMENSION(C%NLON,C%NLAT),               INTENT(IN)  :: lat_gcm
    
    ! Local variables:
    REAL(dp),      DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED)              :: lon_coordinates_of_im_grid_points
    REAL(dp),      DIMENSION(C%NX_EXTENDED,C%NY_EXTENDED)              :: lat_coordinates_of_im_grid_points
    INTEGER,       DIMENSION(C%NLON,C%NLAT)                            :: mask 
    INTEGER                                                            :: i, j
    INTEGER                                                            :: m, n
    INTEGER                                                            :: quadrant          ! The quadrant I, II, III or IV relative to an GCM grid point
    TYPE(triplet)                                                      :: projected_im      ! Projected GCM point on S' 
    TYPE(triplet), DIMENSION(4)                                        :: contribution      ! Nearest projected IM point in quadrant (DIM=I,II,III or IV) in S, relative to the GCM grid point

    ! Projection of the IM coordinates to the GCM coordinates with the inverse oblique stereographic projection:
    ! Output: lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points
    CALL im_projected_lonlat_coordinates(lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
     
    ! A mask is created for the GCM points which are involved (mask ==1) in the mapping:
    ! Output: mask
    CALL make_mask(lon_gcm, lat_gcm, mask)

    ! Opening the 'fast input' file:
    OPEN(UNIT=108, FILE=TRIM(C%input_fast_map_IM_to_GCM))
    
    ! Writing the header of the C%input_fast_map_GCM_to_IM file:
    WRITE(UNIT=108, FMT='(A)') '# Do not remove this header. The quadrant method was used. The format is:'
    WRITE(UNIT=108, FMT='(A)') '#  i  j  N  N(m  n  distance)'
    WRITE(UNIT=108, FMT='(A)') '# with i = the longitudinal GCM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with j = the latitudinal GCM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with N = amount of weighted poits'
    WRITE(UNIT=108, FMT='(A)') '# with m = the x-axis IM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# with n = the y-axis IM grid counter'
    WRITE(UNIT=108, FMT='(A)') '# and distance is the distance between the GCM and IM points'
    WRITE(UNIT=108, FMT='(A, 3(A, F8.4), 2(A, I6))')            '# ', ' alpha_stereographic = ', C%rad2deg * C%alpha_stereographic, ', lambda_M = ', C%rad2deg * C%lambda_M, ', phi_M = ', C%rad2deg * C%phi_M, ', NLON = ', C%NLON, ', NLAT = ', C%NLAT
    WRITE(UNIT=108, FMT='(A, 2(A, I6), 2(A, F12.2), 2(A, I6))') '# ', ' NX = ', C%NX, ', NY = ', C%NY, ', dx = ', C%dx, ', dy = ', C%dy, ', NX_ext = ', C%NX_ext, ', NY_ext = ', C%NY_ext
    WRITE(UNIT=108, FMT='(A)') '# The first number below is the amount of mapped GCM points:'

    ! Number of lines in this 'fast input'  file, equal to the amount mapped IM points:
    WRITE(UNIT=108, FMT='(I10)') COUNT(mask==1)

    ! For each GCM grid point the four nearest projected IM points are determined:
    WRITE(UNIT=*,FMT='(A,I3,A)') '  i =   1 (last i  = ',C%NLON,')'
    DO i = 1, C%NLON
       IF(MOD(i*1.0,10.0) == 0)  WRITE(UNIT=*,FMT='(A,I3)') '  i = ', i
    DO j = 1, C%NLAT
    
      ! Initialize the four nearest quadrant distances to a large value (the width of the IM grid domain):
      contribution(:)%distance = C%dx * C%NX

      IF(mask(i,j) == 1) THEN
       DO m = 1, C%NX_EXTENDED
       DO n = 1, C%NY_EXTENDED
         ! Determine the quadrant in which the projected point lies relative to the considered grid point: 
         ! Output: quadrant
         CALL find_quadrant_around_GCM_grid_point(lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n), &
                                                  lon_gcm(i,j), lat_gcm(i,j), quadrant)

         ! Determine in the curved plane S the distance between the projected IM coordinates relative to the considered GCM grid point:
         projected_im%row_index    = m
         projected_im%column_index = n
         projected_im%distance     = distance_over_curved_surface(lon_gcm(i,j), lat_gcm(i,j), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
         ! In case the projected point coincides with the grid point we put it at the very close distance of 1 centimeter, preventing devision by zero:
         IF(projected_im%distance == 0._dp) projected_im%distance = 0.01_dp

         ! Select the in S projected IM point with the shortest distance to the considered GCM grid point in this quadrant, 
         ! and keep this distance and the IM-grid indices of this IM point in S':
         IF(projected_im%distance < contribution(quadrant)%distance) contribution(quadrant) = projected_im
       END DO
       END DO

       WRITE(UNIT=108,FMT='(3I6, 4(2I6,E23.15))') i, j, 4, &
        contribution(1)%row_index, contribution(1)%column_index, contribution(1)%distance, &
        contribution(2)%row_index, contribution(2)%column_index, contribution(2)%distance, &
        contribution(3)%row_index, contribution(3)%column_index, contribution(3)%distance, &
        contribution(4)%row_index, contribution(4)%column_index, contribution(4)%distance
      END IF
    END DO
    END DO    

    ! Closing the 'fast input' file:
    CLOSE(UNIT=108)
  END SUBROUTINE scan_with_quadrant_method_im_to_gcm



  SUBROUTINE find_quadrant_around_GCM_grid_point(lon_value_projected_point, lat_value_projected_point, &
                                                 lon_value_GCM_grid_point,  lat_value_GCM_grid_point, quadrant)              
    ! Determing the quadrant in which a 'projected point' is situated relative to a GCM grid point.                                           
    !   quadrants:
    !   II  |   I
    !       |     
    !  -----|-----
    !       |     
    !  III  |  IV 
    !
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE             
                              
    ! Input variables:
    REAL(dp), INTENT(IN)  :: lon_value_projected_point
    REAL(dp), INTENT(IN)  :: lat_value_projected_point
    REAL(dp), INTENT(IN)  :: lon_value_GCM_grid_point
    REAL(dp), INTENT(IN)  :: lat_value_GCM_grid_point

    ! Output variable:
    INTEGER,  INTENT(OUT) :: quadrant 

    IF(((lon_value_projected_point - lon_value_GCM_grid_point) <=  180._dp .AND. (lon_value_projected_point - lon_value_GCM_grid_point) > 0._dp) .OR. &
       ((lon_value_projected_point - lon_value_GCM_grid_point) <= -180._dp)) THEN
      IF(lat_value_projected_point > lat_value_GCM_grid_point) THEN
        quadrant = 1
      ELSE
        quadrant = 4
      END IF
    ELSE
      IF(lat_value_projected_point <= lat_value_GCM_grid_point) THEN
        quadrant = 3
      ELSE
        quadrant = 2
      END IF
    END IF
  END SUBROUTINE find_quadrant_around_GCM_grid_point

END MODULE oblimap_scan_contributions_module
