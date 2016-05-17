! File name: oblimap_post_processing_module.f90          
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

MODULE oblimap_post_processing_module

CONTAINS
  SUBROUTINE calculate_average_differences(im_field, gcm_field, initial_gcm_field, mask, field_name, lon_gcm, lat_gcm)
    ! This subroutine calculates the average differences between the initial GCM field and the up and down mapped
    ! field, aslo the standard deviations are calculated, and their distribution written to a file for plotting.  
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NY,  C%NX  ), INTENT(IN) :: im_field
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN) :: gcm_field
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN) :: initial_gcm_field
    INTEGER,  DIMENSION(C%NLON,C%NLAT), INTENT(IN) :: mask
    CHARACTER(LEN=*),                   INTENT(IN) :: field_name
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN) :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN) :: lat_gcm
        
    ! Local variables:
    REAL(dp), DIMENSION(C%NLON, C%NLAT)            :: gcm_field_differences
    REAL(dp)                                       :: lowest_gcm_field_value 
    REAL(dp)                                       :: highest_gcm_field_value
    REAL(dp)                                       :: average_difference
    REAL(dp)                                       :: lowest_difference
    REAL(dp)                                       :: highest_difference
    REAL(dp)                                       :: gcm_range
    REAL(dp)                                       :: average_mapped    
    REAL(dp)                                       :: standard_deviation     
    CHARACTER(LEN=255)                             :: distribution_file_name
    REAL(dp)                                       :: delta    
    REAL(dp)                                       :: a    
    INTEGER                                        :: N
    INTEGER                                        :: counter     
    CHARACTER(LEN=128)                             :: data_label
    CHARACTER(LEN=256)                             :: print_format
    REAL(dp)                                       :: gcm_area_by_triangles            ! gcm area determined by the triangle method
    REAL(dp)                                       :: gcm_area_integrated_field_value_by_triangles  
    REAL(dp)                                       :: area_weighted_field_average_gcm_by_triangles
    REAL(dp)                                       :: area_weighted_field_average_im

    ! Output: gcm_area_by_triangles, gcm_area_integrated_field_value_by_triangles
    CALL calculate_integrated_quantity_gcm_area_by_triangles(  lon_gcm, lat_gcm, initial_gcm_field, mask, gcm_area_by_triangles,   gcm_area_integrated_field_value_by_triangles)
    area_weighted_field_average_gcm_by_triangles   = gcm_area_integrated_field_value_by_triangles   / gcm_area_by_triangles      ! area weighted values gcm for mask==1
    area_weighted_field_average_im                 = SUM(im_field) / (C%NX * C%NY)                                               ! area weighted values im

    gcm_field_differences   = gcm_field - initial_gcm_field
    lowest_gcm_field_value  = MINVAL(initial_gcm_field, mask == 1)
    highest_gcm_field_value = MAXVAL(initial_gcm_field, mask == 1)
    gcm_range               = (highest_gcm_field_value - lowest_gcm_field_value)
    average_mapped          = SUM(initial_gcm_field, mask == 1) / COUNT(mask == 1)
    average_difference      = SUM(ABS(gcm_field_differences), mask == 1) / COUNT(mask == 1)
    lowest_difference       = MINVAL(gcm_field_differences, mask == 1)
    highest_difference      = MAXVAL(gcm_field_differences, mask == 1)
    standard_deviation      = (SUM((gcm_field_differences - average_difference)**2, mask == 1) / COUNT(mask == 1))**0.5_dp
    
    WRITE(UNIT=8, FMT='(2(A, F12.4), 8(A, F12.6), A, A10)') &
     'M = (lon='                , C%rad2deg * C%lambda_M                       , &
     ',lat='                    , C%rad2deg * C%phi_M                          , &
     '), Average (at IM) = '    , area_weighted_field_average_im               , &
     ', Average difference = '  , average_difference                           , &
     ' = '                      , 100._dp * average_difference / average_mapped, &
     '%, Lowest difference = '  , lowest_difference                            , &
     ', Hightest difference = ' , highest_difference                           , &
     ', standard deviation = '  , standard_deviation                           , &
     ', Within '                , 2._dp * standard_deviation                   , &
     ', being ', 100._dp * COUNT(ABS(gcm_field_differences) < 2._dp * standard_deviation .AND. mask == 1)/COUNT(mask == 1), &
     '% confidence, '           , field_name
     
    print_format  = '(2(A, F9.1), 3(A, F9.1), 2(A, F12.2), 2(A, F12.2), 2(A, F12.3), 2(A, ES11.2), 2(A, F9.2), (A, F10.2), (A, F14.5), (A, F7.1), 3(A, E11.2), 2(A, I5), (A, F4.1), 3(A, F8.1), (A,I10), 2A, A10)'
    IF(field_name == 'Ts') &    
     print_format = '(2(A, F9.1), 3(A, F9.1), 2(A, F12.2), 2(A, F12.2), 2(A, F12.3), 2(A, ES11.2), 2(A, F9.2), (A, F10.2), (A, F14.5), (A, F7.1), 3(A, E11.2), 2(A, I5), (A, F4.1), 3(A, F8.1), (A,I10), 2A, A10)'
    IF(field_name == 'MB_surface') &    
     print_format = '(2(A, F9.2), 3(A, F9.2), 2(A, F12.3), 2(A, F12.2), 2(A, F12.4), 2(A, ES11.2), 2(A, F9.3), (A, F10.2), (A, F14.5), (A, F7.1), 3(A, E11.2), 2(A, I5), (A, F4.1), 3(A, F8.1), (A,I10), 2A, A10)'
    IF(field_name == 'Hs') &    
     print_format = '(2(A, F9.0), 3(A, F9.0), 2(A, F12.1), 2(A, F12.2), 2(A, F12.1), 2(A, ES11.2), 2(A, F9.1), (A, F10.2), (A, F14.5), (A, F7.1), 3(A, E11.2), 2(A, I5), (A, F4.1), 3(A, F8.1), (A,I10), 2A, A10)'
    ! Write in latex table format:
    WRITE(UNIT=81, FMT=print_format) &
     ' & ', lowest_gcm_field_value                                                                                    , & ! lowest GCM value of the considered quantity
     ' & ', highest_gcm_field_value                                                                                   , & ! highest GCM value of the considered quantity
     ' &R', average_mapped                                                                                            , & ! point average, this average is not always area-weighted
     ' &C', area_weighted_field_average_gcm_by_triangles                                                              , & ! area-weigted average over GCM grid by the triangle method
     ' &:', area_weighted_field_average_im                                                                            , & ! area-weigted average over IM grid
     ' &R',           ABS(average_mapped                                 - area_weighted_field_average_im)            , & ! difference point        average over GCM grid                           with area-weigted average over IM grid
     ' &C',           ABS(area_weighted_field_average_gcm_by_triangles   - area_weighted_field_average_im)            , & ! difference area-weigted average over GCM grid by the triangle method    with area-weigted average over IM grid
     ' &R', 100._dp * ABS(area_weighted_field_average_im - average_mapped                                ) / gcm_range, & ! rrd area-weigted averages IM - point        average over GCM grid                           in percent
     ' &C', 100._dp * ABS(area_weighted_field_average_im - area_weighted_field_average_gcm_by_triangles  ) / gcm_range, & ! rrd area-weigted averages IM - area-weigted average over GCM grid by the triangle method    in percent
     ' &:', 100._dp * ABS(1._dp - area_weighted_field_average_im / average_mapped                                )    , & ! fraction area-weigted averages IM / point        average over GCM grid                           in percent
     ' &:', 100._dp * ABS(1._dp - area_weighted_field_average_im / area_weighted_field_average_gcm_by_triangles  )    , & ! fraction area-weigted averages IM / area-weigted average over GCM grid by the triangle method    in percent
     ' &R', ABS(average_mapped                                 - area_weighted_field_average_im) * C%NX*C%NY*C%dx*C%dy, & ! difference point        average over GCM grid                           with area-weigted average over IM grid
     ' &C', ABS(area_weighted_field_average_gcm_by_triangles   - area_weighted_field_average_im) * C%NX*C%NY*C%dx*C%dy, & ! difference area-weigted average over GCM grid by the triangle method    with area-weigted average over IM grid
     ' & ', average_difference                                                                                        , & ! point average of the difference, this average is not area-weighted
     ' & ', 2._dp * standard_deviation                                                                                , & ! two sigma, sima is the standard deviation
     ' & ', 100._dp * average_difference / gcm_range                                                                  , & ! range relative deviation (RRD)
     ' &:', 100._dp * SUM(ABS(gcm_field_differences / initial_gcm_field), mask == 1 .AND. initial_gcm_field /= 0._dp) / COUNT(mask == 1 .AND. initial_gcm_field /= 0._dp), &
     ' &:', 100._dp * COUNT(ABS(gcm_field_differences) < 2._dp * standard_deviation .AND. mask == 1)/COUNT(mask == 1) , & ! percentage of points which lie within two sigma
     ' &:', 100._dp * average_difference / average_mapped                                                             , & ! difference in percent relative to the average
     ' &:', lowest_difference                                                                                         , &
     ' &:', highest_difference                                                                                        , &
     ' &:', C%NX                                                                                                      , &
     ' &:', C%NY                                                                                                      , &
     ' &:', C%dx / 1000._dp                                                                                           , &
     ' &:', C%rad2deg * C%alpha_stereographic                                                                         , &
     ' &:', C%rad2deg * C%lambda_M                                                                                    , &
     ' &:', C%rad2deg * C%phi_M                                                                                       , &
     ' &:', COUNT(mask == 1), ' \\\\'                                                                                 , &
     '   ', field_name

    ! Writing the basic file for the distribution plots:
    
    distribution_file_name = 'distribution_'//TRIM(field_name)//'.dat'
    OPEN(UNIT=9000, FILE=distribution_file_name)

    N = 300
    ! Showing all points, also the far outlayers:
    delta = (highest_difference - lowest_difference) / REAL(N, dp)
    a = lowest_difference
    ! Showing all points within two sigma:
    !delta = 2._dp * standard_deviation / REAL(N, dp)
    !a = - standard_deviation
    DO counter = 1, N
     IF(ABS(a) > 2._dp * standard_deviation) THEN
      data_label = '  ! data outside sigma x 2'
     ELSE IF(ABS(a) > standard_deviation) THEN
      data_label = '  ! data inside sigma x 2'
     ELSE
      data_label = '  ! data inside sigma x 1 '
     END IF
     WRITE(UNIT=9000, FMT='(F12.6, I10, A)') a, COUNT(gcm_field_differences >= a .AND. gcm_field_differences < a + delta .AND. mask == 1), data_label
     a = a + delta
    END DO
    
    CLOSE(UNIT=9000)
  END SUBROUTINE calculate_average_differences



  ! Under development:
  ! The GCM grid edges should be taken more proper an consistent, although these cases rarely occur (but can be improved).
  SUBROUTINE calculate_integrated_quantity_gcm_area_by_triangles(lon_gcm, lat_gcm, gcm_field, mask, gcm_area, gcm_area_integrated_value)
    ! This subroutine calculates the area integrated value of a quantity at the gcm points which
    ! participate in the mapping, i.e. with mask == 1.
    !
    ! Determine the four corner points A, B, C, and D of the grid box around lon_gcm(i,j), lat_gcm(i,j) 
    ! The sketch below shows the location of A = lambda_a, phi_a, B = lambda_b, phi_b, etc.
    ! You can see that for instance that D is the center point of (lon_gcm(i,j), lat_gcm(i,j)),
    ! (lon_gcm(i+1,j), lat_gcm(i+1,j)), (lon_gcm(i,j+1), lat_gcm(i,j+1)) and (lon_gcm(i+1,j+1), lat_gcm(i+1,j+1))
    !
    !                                      (i+1, j+1)
    !                                     x
    !
    !                    (i,j+1)
    !                           x     
    !                                D        
    !                                 .         x (i+1,j)
    !
    !             x                           
    !                     C.         
    !                              x        .
    !                             (i,j)     B       x
    !                 x                                
    !                          .                        
    !                          A        x
    !                                         
    !                       x
    !
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: gcm_field
    INTEGER,  DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: mask

    ! Output variable:
    REAL(dp),                           INTENT(OUT) :: gcm_area
    REAL(dp),                           INTENT(OUT) :: gcm_area_integrated_value
        
    ! Local variables:
    INTEGER                                         :: i, j
    INTEGER                                         :: im1  ! = i - 1
    INTEGER                                         :: ip1  ! = i + 1
    INTEGER                                         :: jm1  ! = j - 1
    INTEGER                                         :: jp1  ! = j + 1

    REAL(dp)                                        :: lambda_a, lambda_b, lambda_c, lambda_d
    REAL(dp)                                        :: phi_a, phi_b, phi_c, phi_d
    REAL(dp)                                        :: area
    
    gcm_area                  = 0._dp
    gcm_area_integrated_value = 0._dp
    DO i = 1, C%NLON 
    DO j = 1, C%NLAT 
      IF(mask(i,j) == 1) THEN
       ! Determine if the grid of the initial GCM file represents the whole earth or only a region:
       IF (lon_gcm(i,j) + lon_gcm(C%NLON,j) + (lon_gcm(C%NLON,j) - lon_gcm(C%NLON - 1,j)) == 360) THEN
        ! GCM grid is whole earth -> Last longitude for i=C%NLON is followed by i=1 again.
        IF (i == 1) THEN
          im1 = C%NLON
        ELSE
          im1 = i - 1
        ENDIF
        IF (i == C%NLON) THEN
          ip1 = 1
        ELSE
          ip1 = i + 1
        ENDIF
       ELSE
        IF (i == 1) THEN
          im1 = i
        ELSE
          im1 = i - 1
        ENDIF
        IF (i == C%NLON) THEN
          ip1 = i
        ELSE
          ip1 = i + 1
        ENDIF
       ENDIF

       IF (j == 1) THEN
        jm1 = j
       ELSE
        jm1 = j - 1
       ENDIF

       IF (j == C%NLAT) THEN
        jp1 = j
       ELSE
        jp1 = j + 1
       ENDIF

       ! Determine the four corner points of the grid box surrounding point (i,j):
       CALL find_center_point(lon_gcm(im1,jm1), lat_gcm(im1,jm1),lon_gcm(i  ,jm1), lat_gcm(i  ,jm1),lon_gcm(im1,j  ), lat_gcm(im1,j  ),lon_gcm(i  ,j  ), lat_gcm(i  ,j  ),lambda_a, phi_a)
       CALL find_center_point(lon_gcm(i  ,jm1), lat_gcm(i  ,jm1),lon_gcm(ip1,jm1), lat_gcm(ip1,jm1),lon_gcm(i  ,j  ), lat_gcm(i  ,j  ),lon_gcm(ip1,j  ), lat_gcm(ip1,j  ),lambda_b, phi_b)
       CALL find_center_point(lon_gcm(im1,j  ), lat_gcm(im1,j  ),lon_gcm(i  ,j  ), lat_gcm(i  ,j  ),lon_gcm(im1,jp1), lat_gcm(im1,jp1),lon_gcm(i  ,jp1), lat_gcm(i  ,jp1),lambda_c, phi_c)
       CALL find_center_point(lon_gcm(i  ,j  ), lat_gcm(i  ,j  ),lon_gcm(ip1,j  ), lat_gcm(ip1,j  ),lon_gcm(i  ,jp1), lat_gcm(i  ,jp1),lon_gcm(ip1,jp1), lat_gcm(ip1,jp1),lambda_d, phi_d)

       ! Determine area between four corner points  
       CALL area_grid_box(lambda_a, phi_a, lambda_b, phi_b, lambda_c, phi_c, lambda_d, phi_d, area)
       gcm_area                  = gcm_area                  + area
       gcm_area_integrated_value = gcm_area_integrated_value + area * gcm_field(i,j)
      END IF
    END DO 
    END DO

    ! Adhoc calulation of the total RACMO area:
    !CALL area_grid_box(303.517_dp,54.3987_dp, 344.257_dp,53.8027_dp, 244.768_dp,78.7983_dp, 42.3487_dp,77.1189_dp, area)
    !WRITE(6,*) 'The total RACMO area =', area
  END SUBROUTINE calculate_integrated_quantity_gcm_area_by_triangles



  ! Under development:
  ! Actually the center point should be determined by finding the intersection point of the
  ! two diagonals over the great circles (currently we took the local surface flat)
  SUBROUTINE find_center_point(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4, lonc, latc)
    ! This subroutine determines the central point C between the four points 1, 2, 3, and 4:
    !                 
    !                 x 4
    !     3  x
    !               C
    !                       x 2
    !           1  x
    !                    
    ! lon1,lat1 = lon_gcm(i  ,j  ),lat_gcm(i  ,j  )                                            
    ! lon2,lat2 = lon_gcm(i+1,j  ),lat_gcm(i+1,j  )                                            
    ! lon3,lat3 = lon_gcm(i  ,j+1),lat_gcm(i  ,j+1)                                            
    ! lon4,lat4 = lon_gcm(i+1,j+1),lat_gcm(i+1,j+1) 
    ! lonc,latc = coordinates of the central point C                                         
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lon1, lat1
    REAL(dp), INTENT(IN)  :: lon2, lat2
    REAL(dp), INTENT(IN)  :: lon3, lat3
    REAL(dp), INTENT(IN)  :: lon4, lat4

    ! Output variables:
    REAL(dp),  INTENT(OUT) :: lonc, latc

    ! Local variables:
    ! which follow from the point of intersection between the two diagonals 
    ! d1 = diagonal between point 1 and 4 
    ! d2 = diagonal between point 2 and 3
    ! both represented by a parameter equation
    REAL(dp)  :: r12  ! parameter that describes certain relation between point 1 and 2
    REAL(dp)  :: r34  ! parameter that describes certain relation between point 3 and 4
    REAL(dp)  :: lon2_n
    REAL(dp)  :: lon4_n
   
    ! Central point C(lonc, latc) is given by following formulas:
    if (lon1 == lon2 .AND. lon3 == lon4) then
     lonc = lon1 + 0.5_dp * (lon3 - lon1)
     latc = lat1 + 0.5_dp * (lat3 - lat1)
    else
     if (lon1 > lon2) then
      lon2_n = lon2 + 360_dp
     else
      lon2_n = lon2 
     endif
     if (lon3 > lon4) then
      lon4_n = lon4 + 360_dp
     else
      lon4_n = lon4 
     endif
     r12 = (lon2_n - lon1) + (lat2 - lat1)
     r34 = (lon4_n - lon3) + (lat4 - lat3)

     if (lonc > 360_dp) then
      lonc = lonc - 360_dp
     else
      lonc = lon1 + (lon4_n - lon1) * (r12 / (r12 + r34))    
     endif
     latc = lat1 + (lat4 - lat1) * (r12 / (r12 + r34))    
    endif
  END SUBROUTINE find_center_point



 SUBROUTINE area_grid_box(lambda_a, phi_a, lambda_b, phi_b, lambda_c, phi_c, lambda_d, phi_d, area)
    ! This subroutine determines the area of grid box X which lies between the four circled (o) edge points A, B, C, and D:
    !
    !                 o D
    !     C  o
    !               X
    !                       o B
    !           A  o
    !                    
    ! area = the area of grid box X                                           
    !
    ! To calculate the area of grid box X, we determine the areas of the two triangles: area = area ACD + area ABD
    ! The surface of one triangle is calulated by trigoniometry on a spherical surface: Finding the area of a spherical triangle
    ! See: http://mathforum.org/library/drmath/view/65316.html
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_scan_contributions_module, ONLY: distance_over_curved_surface
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda_a, phi_a
    REAL(dp), INTENT(IN)  :: lambda_b, phi_b
    REAL(dp), INTENT(IN)  :: lambda_c, phi_c
    REAL(dp), INTENT(IN)  :: lambda_d, phi_d

    ! Output variable:
    REAL(dp), INTENT(OUT) :: area

    ! Local variables:
    REAL(dp)              :: d_AC, d_CD, d_AD, d_AB, d_BD
    REAL(dp)              :: s_upper, s_lower
    REAL(dp)              :: E_upper, E_lower
    REAL(dp)              :: area_upper, area_lower

    ! First determine all the distances over the sphere of all the sides: 
    d_AC = distance_over_curved_surface(lambda_a, phi_a, lambda_c, phi_c) / C%R_earth
    d_CD = distance_over_curved_surface(lambda_c, phi_c, lambda_d, phi_d) / C%R_earth
    d_AD = distance_over_curved_surface(lambda_a, phi_a, lambda_d, phi_d) / C%R_earth
    d_AB = distance_over_curved_surface(lambda_a, phi_a, lambda_b, phi_b) / C%R_earth
    d_BD = distance_over_curved_surface(lambda_b, phi_b, lambda_d, phi_d) / C%R_earth

    ! Determine the area of the upper triangle { ACD } 
    s_upper    = ( d_AC + d_CD + d_AD ) / 2._dp
    E_upper    = 4._dp * ATAN(SQRT(TAN(s_upper / 2._dp) * TAN((s_upper - d_AC) / 2._dp) * TAN((s_upper - d_CD) / 2._dp) * TAN((s_upper - d_AD) / 2._dp)))
    area_upper = C%R_earth**2 * E_upper

    ! Determine the area of the upper triangle { ABD } 
    s_lower    = ( d_AB + d_AD + d_BD ) / 2._dp
    E_lower    = 4._dp * ATAN(SQRT(TAN(s_lower / 2._dp) * TAN((s_lower - d_AB) / 2._dp) * TAN((s_lower - d_AD) / 2._dp) * TAN((s_lower - d_BD) / 2._dp)))
    area_lower = C%R_earth**2 * E_lower

    ! Determine the area of grid box X: 
    area       = area_upper + area_lower
  END SUBROUTINE area_grid_box
      
END MODULE oblimap_post_processing_module

