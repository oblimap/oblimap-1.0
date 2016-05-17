! File name: oblimap_read_and_write_module.f90          
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

MODULE oblimap_read_and_write_module

CONTAINS
  SUBROUTINE read_data_gcm(file_gcm, field_name_1_gcm,  lon_gcm, lat_gcm, field_1_gcm, &
                                     field_name_2_gcm,                    field_2_gcm, &
                                     field_name_3_gcm,                    field_3_gcm, &
                                     field_name_4_gcm,                    field_4_gcm, &
                                     field_name_5_gcm,                    field_5_gcm, &
                                     field_name_6_gcm,                    field_6_gcm, &
                                     field_name_7_gcm,                    field_7_gcm, &
                                     field_name_8_gcm,                    field_8_gcm, &
                                     field_name_9_gcm,                    field_9_gcm)
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                   INTENT(IN)            :: file_gcm         ! GCM input file
    CHARACTER(LEN=*),                   INTENT(IN)            :: field_name_1_gcm ! name of a 1st GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_2_gcm ! name of a 2nd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_3_gcm ! name of a 3rd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_4_gcm ! name of a 4th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_5_gcm ! name of a 5th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_6_gcm ! name of a 6th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_7_gcm ! name of a 7th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_8_gcm ! name of a 8th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_9_gcm ! name of a 9th GCM variable

    ! Output variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: lon_gcm          ! longitude coordinates (in degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: lat_gcm          ! latitude coordinates  (in degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: field_1_gcm      ! values of a 1st GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_2_gcm      ! values of a 2nd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_3_gcm      ! values of a 3rd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_4_gcm      ! values of a 4th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_5_gcm      ! values of a 5th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_6_gcm      ! values of a 6th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_7_gcm      ! values of a 7th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_8_gcm      ! values of a 8th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_9_gcm      ! values of a 9th GCM field 

    ! Local variables:
    REAL(dp), DIMENSION(C%NLON)                               :: lon_gcm_1D       ! 1D longitude coordinates (in degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLAT)                               :: lat_gcm_1D       ! 1D latitude coordinates  (in degrees) of GCM grid
    INTEGER                                                   :: ncid            
    INTEGER                                                   :: lon_id_1D       
    INTEGER                                                   :: lat_id_1D       
    INTEGER                                                   :: var_id          
    INTEGER                                                   :: i, j            
                                                                                 
    ! Read the variables:
    CALL handle_error(nf90_open(file_gcm, nf90_nowrite, ncid))
    CALL handle_error(nf90_inq_varid(ncid, C%gcm_x_axis_name, lon_id_1D))
    CALL handle_error(nf90_inq_varid(ncid, C%gcm_y_axis_name, lat_id_1D))
    CALL handle_error(nf90_inq_varid(ncid, field_name_1_gcm, var_id))
    CALL handle_error(nf90_get_var(ncid, lon_id_1D, lon_gcm_1D, start=(/ C%starting_recordnr/)))
    CALL handle_error(nf90_get_var(ncid, lat_id_1D, lat_gcm_1D, start=(/ C%starting_recordnr/)))
    CALL handle_error(nf90_get_var(ncid, var_id, field_1_gcm, start=(/1, 1, C%starting_recordnr/)))
          
    IF(PRESENT(field_2_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_2_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_2_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_3_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_3_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_3_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_4_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_4_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_4_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_5_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_5_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_5_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_6_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_6_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_6_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_7_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_7_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_7_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_8_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_8_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_8_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_9_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_9_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_9_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    
    ! Creating from the one dimensional lon lat variables, the two dimensional ones:
    DO i = 1, C%NLON
     lon_gcm(i,:) = lon_gcm_1D(i)
    END DO
    DO j = 1, C%NLAT
     lat_gcm(:,j) = lat_gcm_1D(j)
    END DO
          
    ! Close the file:
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE read_data_gcm



  SUBROUTINE read_data_gcm_with_2D_lonlat(file_gcm, field_name_1_gcm, lon_gcm, lat_gcm, field_1_gcm, &
                                                    field_name_2_gcm,                   field_2_gcm, &
                                                    field_name_3_gcm,                   field_3_gcm, &
                                                    field_name_4_gcm,                   field_4_gcm, &
                                                    field_name_5_gcm,                   field_5_gcm, &
                                                    field_name_6_gcm,                   field_6_gcm, &
                                                    field_name_7_gcm,                   field_7_gcm, &
                                                    field_name_8_gcm,                   field_8_gcm, &
                                                    field_name_9_gcm,                   field_9_gcm)
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                   INTENT(IN)            :: file_gcm         ! GCM input file
    CHARACTER(LEN=*),                   INTENT(IN)            :: field_name_1_gcm ! name of a 1st GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_2_gcm ! name of a 2nd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_3_gcm ! name of a 3rd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_4_gcm ! name of a 4th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_5_gcm ! name of a 5th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_6_gcm ! name of a 6th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_7_gcm ! name of a 7th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_8_gcm ! name of a 8th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN),  OPTIONAL :: field_name_9_gcm ! name of a 9th GCM variable

    ! Output variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: lon_gcm          ! longitude coordinates (in degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: lat_gcm          ! latitude coordinates  (in degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT)           :: field_1_gcm      ! values of a 1st GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_2_gcm      ! values of a 2nd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_3_gcm      ! values of a 3rd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_4_gcm      ! values of a 4th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_5_gcm      ! values of a 5th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_6_gcm      ! values of a 6th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_7_gcm      ! values of a 7th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_8_gcm      ! values of a 8th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT), OPTIONAL :: field_9_gcm      ! values of a 9th GCM field 

    ! Local variables:
    INTEGER                                                   :: ncid
    INTEGER                                                   :: lon_id
    INTEGER                                                   :: lat_id
    INTEGER                                                   :: var_id

    ! Read the variables:
    CALL handle_error(nf90_open(file_gcm, nf90_nowrite, ncid))
    CALL handle_error(nf90_inq_varid(ncid, C%gcm_x_axis_name, lon_id))
    CALL handle_error(nf90_inq_varid(ncid, C%gcm_y_axis_name, lat_id))
    CALL handle_error(nf90_inq_varid(ncid, field_name_1_gcm, var_id))
    CALL handle_error(nf90_get_var(ncid, lon_id, lon_gcm, start=(/1, 1, C%starting_recordnr/)))
    CALL handle_error(nf90_get_var(ncid, lat_id, lat_gcm, start=(/1, 1, C%starting_recordnr/)))
    CALL handle_error(nf90_get_var(ncid, var_id, field_1_gcm, start=(/1, 1, C%starting_recordnr/)))
          
    IF(PRESENT(field_2_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_2_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_2_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_3_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_3_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_3_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_4_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_4_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_4_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_5_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_5_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_5_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_6_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_6_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_6_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_7_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_7_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_7_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_8_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_8_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_8_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_9_gcm)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_9_gcm,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_9_gcm, start=(/1, 1, C%starting_recordnr/)))
    END IF
    
    ! Close the file:
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE read_data_gcm_with_2D_lonlat


  
  SUBROUTINE write_data_for_im(file_im, field_name_1_im, field_1_im, &
                                        field_name_2_im, field_2_im, &
                                        field_name_3_im, field_3_im, &
                                        field_name_4_im, field_4_im, &
                                        field_name_5_im, field_5_im, &
                                        field_name_6_im, field_6_im, &
                                        field_name_7_im, field_7_im, &
                                        field_name_8_im, field_8_im, &
                                        field_name_9_im, field_9_im)
    ! This routine writes projected + interpolated IM field to a netcdf file.
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),               INTENT(IN)            :: file_im
    CHARACTER(LEN=*),               INTENT(IN)            :: field_name_1_im ! name of a 1st IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_2_im ! name of a 2nd IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_3_im ! name of a 3rd IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_4_im ! name of a 4th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_5_im ! name of a 5th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_6_im ! name of a 6th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_7_im ! name of a 7th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_8_im ! name of a 8th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_9_im ! name of a 9th IM variable

    ! Output variable:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(IN)            :: field_1_im      ! values of a 1st IM field
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_2_im      ! values of a 2nd IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_3_im      ! values of a 3rd IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_4_im      ! values of a 4th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_5_im      ! values of a 5th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_6_im      ! values of a 6th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_7_im      ! values of a 7th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_8_im      ! values of a 8th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_9_im      ! values of a 9th IM field 
    
    ! Local variables:
    REAL(dp), DIMENSION(C%NX)                             :: x_im
    REAL(dp), DIMENSION(C%NY)                             :: y_im
    INTEGER                                               :: m, n
    INTEGER                                               :: x_dim, y_dim
    INTEGER                                               :: ncid
    INTEGER                                               :: x_id, y_id
    INTEGER                                               :: var_1_id
    INTEGER                                               :: var_2_id
    INTEGER                                               :: var_3_id
    INTEGER                                               :: var_4_id
    INTEGER                                               :: var_5_id
    INTEGER                                               :: var_6_id
    INTEGER                                               :: var_7_id
    INTEGER                                               :: var_8_id
    INTEGER                                               :: var_9_id
    
    ! Determine the x-y coordinates of the IM gridpoints on projection plane P
    ! In total there are C%NX x C%NY gridpoints and the gridcell distances are C%dx and C%dy (in meters).
    ! The center gridpoint is always (0,0) and represents the longitude-latitude coordinates
    ! (lamda_M,phi_M). For instance if the region of interest is Antarctica, then
    ! (lamda_M,phi_M) = (0,-90) representing the Southpole.

    DO m = 1, C%NX
      x_im(m) = C%dx * (m - ((C%NX+1) / 2))
    END DO
    
    DO n = 1, C%NY
      y_im(n) = C%dy * (n - ((C%NY+1) / 2))
    END DO
    
    ! Write variables to file:
    CALL handle_error(nf90_create(file_im,nf90_clobber,ncid))
    ! Define dimensions 
    CALL handle_error(nf90_def_dim(ncid, 'X', C%NX, x_dim))
    CALL handle_error(nf90_def_dim(ncid, 'Y', C%NY, y_dim))
    ! Define variables
    CALL handle_error(nf90_def_var(ncid, 'X', nf90_int, (/x_dim/), x_id))
    CALL handle_error(nf90_def_var(ncid, 'Y', nf90_int, (/y_dim/), y_id))
    CALL handle_error(nf90_def_var(ncid, field_name_1_im, nf90_float, (/x_dim, y_dim/), var_1_id))
    IF(PRESENT(field_2_im)) CALL handle_error(nf90_def_var(ncid, field_name_2_im , nf90_float, (/x_dim, y_dim/), var_2_id))
    IF(PRESENT(field_3_im)) CALL handle_error(nf90_def_var(ncid, field_name_3_im , nf90_float, (/x_dim, y_dim/), var_3_id))
    IF(PRESENT(field_4_im)) CALL handle_error(nf90_def_var(ncid, field_name_4_im , nf90_float, (/x_dim, y_dim/), var_4_id))
    IF(PRESENT(field_5_im)) CALL handle_error(nf90_def_var(ncid, field_name_5_im , nf90_float, (/x_dim, y_dim/), var_5_id))
    IF(PRESENT(field_6_im)) CALL handle_error(nf90_def_var(ncid, field_name_6_im , nf90_float, (/x_dim, y_dim/), var_6_id))
    IF(PRESENT(field_7_im)) CALL handle_error(nf90_def_var(ncid, field_name_7_im , nf90_float, (/x_dim, y_dim/), var_7_id))
    IF(PRESENT(field_8_im)) CALL handle_error(nf90_def_var(ncid, field_name_8_im , nf90_float, (/x_dim, y_dim/), var_8_id))
    IF(PRESENT(field_9_im)) CALL handle_error(nf90_def_var(ncid, field_name_9_im , nf90_float, (/x_dim, y_dim/), var_9_id))
    CALL handle_error(nf90_enddef(ncid))
    CALL handle_error(nf90_put_var(ncid, x_id,   x_im,     start=(/      1/)))
    CALL handle_error(nf90_put_var(ncid, y_id,   y_im,     start=(/      1/)))
    CALL handle_error(nf90_put_var(ncid, var_1_id, field_1_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_2_im)) CALL handle_error(nf90_put_var(ncid, var_2_id, field_2_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_3_im)) CALL handle_error(nf90_put_var(ncid, var_3_id, field_3_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_4_im)) CALL handle_error(nf90_put_var(ncid, var_4_id, field_4_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_5_im)) CALL handle_error(nf90_put_var(ncid, var_5_id, field_5_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_6_im)) CALL handle_error(nf90_put_var(ncid, var_6_id, field_6_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_7_im)) CALL handle_error(nf90_put_var(ncid, var_7_id, field_7_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_8_im)) CALL handle_error(nf90_put_var(ncid, var_8_id, field_8_im, start=(/1, 1, 1/)))
    IF(PRESENT(field_9_im)) CALL handle_error(nf90_put_var(ncid, var_9_id, field_9_im, start=(/1, 1, 1/)))

    ! close the file 
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE write_data_for_im



  SUBROUTINE read_data_im(file_im, field_name_1_im, field_1_im, &
                                   field_name_2_im, field_2_im, &
                                   field_name_3_im, field_3_im, &
                                   field_name_4_im, field_4_im, &
                                   field_name_5_im, field_5_im, &
                                   field_name_6_im, field_6_im, &
                                   field_name_7_im, field_7_im, &
                                   field_name_8_im, field_8_im, &
                                   field_name_9_im, field_9_im)
    ! This routine reads the IM field(s)
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),               INTENT(IN)            :: file_im 
    CHARACTER(LEN=*),               INTENT(IN)            :: field_name_1_im ! name of a 1st IM variable 
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_2_im ! name of a 2nd IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_3_im ! name of a 3rd IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_4_im ! name of a 4th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_5_im ! name of a 5th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_6_im ! name of a 6th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_7_im ! name of a 7th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_8_im ! name of a 8th IM variable
    CHARACTER(LEN=*),               INTENT(IN),  OPTIONAL :: field_name_9_im ! name of a 9th IM variable
                                   
    ! Output variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT)           :: field_1_im      ! values of a 1st IM field
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_2_im      ! values of a 2nd IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_3_im      ! values of a 3rd IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_4_im      ! values of a 4th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_5_im      ! values of a 5th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_6_im      ! values of a 6th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_7_im      ! values of a 7th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_8_im      ! values of a 8th IM field 
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT), OPTIONAL :: field_9_im      ! values of a 9th IM field 

    ! Local variables:
    INTEGER                                     :: ncid
    INTEGER                                     :: var_id
  
    ! Read the variables:
    CALL handle_error(nf90_open(file_im,nf90_nowrite,ncid))
    CALL handle_error(nf90_inq_varid(ncid, field_name_1_im, var_id))
    CALL handle_error(nf90_get_var(ncid, var_id, field_1_im, start=(/1, 1, 1/)))

    IF(PRESENT(field_2_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_2_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_2_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_3_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_3_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_3_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_4_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_4_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_4_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_5_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_5_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_5_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_6_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_6_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_6_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_7_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_7_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_7_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_8_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_8_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_8_im, start=(/1, 1, C%starting_recordnr/)))
    END IF
    IF(PRESENT(field_9_im)) THEN
     CALL handle_error(nf90_inq_varid(ncid, field_name_9_im,    var_id))
     CALL handle_error(nf90_get_var(ncid, var_id, field_9_im, start=(/1, 1, C%starting_recordnr/)))
    END IF

    ! close the file
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE read_data_im



  SUBROUTINE write_data_for_gcm(file_gcm, lon_gcm, lat_gcm, field_name_1_gcm, field_1_gcm, &
                                                            field_name_2_gcm, field_2_gcm, &
                                                            field_name_3_gcm, field_3_gcm, &
                                                            field_name_4_gcm, field_4_gcm, &
                                                            field_name_5_gcm, field_5_gcm, &
                                                            field_name_6_gcm, field_6_gcm, &
                                                            field_name_7_gcm, field_7_gcm, &
                                                            field_name_8_gcm, field_8_gcm, &
                                                            field_name_9_gcm, field_9_gcm)
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                   INTENT(IN)           :: file_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lat_gcm
    CHARACTER(LEN=*),                   INTENT(IN)           :: field_name_1_gcm ! name of a 1st GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_2_gcm ! name of a 2nd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_3_gcm ! name of a 3rd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_4_gcm ! name of a 4th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_5_gcm ! name of a 5th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_6_gcm ! name of a 6th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_7_gcm ! name of a 7th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_8_gcm ! name of a 8th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_9_gcm ! name of a 9th GCM variable
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: field_1_gcm      ! values of a 1st GCM field
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_2_gcm      ! values of a 2nd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_3_gcm      ! values of a 3rd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_4_gcm      ! values of a 4th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_5_gcm      ! values of a 5th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_6_gcm      ! values of a 6th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_7_gcm      ! values of a 7th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_8_gcm      ! values of a 8th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_9_gcm      ! values of a 9th GCM field 
    
    ! Local variables:
    REAL(dp), DIMENSION(C%NLON)                              :: lon_gcm_1D 
    REAL(dp), DIMENSION(C%NLAT)                              :: lat_gcm_1D 
    INTEGER                                                  :: ncid
    INTEGER                                                  :: lon_dim
    INTEGER                                                  :: lat_dim
    INTEGER                                                  :: lon_id_1D
    INTEGER                                                  :: lat_id_1D
    INTEGER                                                  :: var_1_id
    INTEGER                                                  :: var_2_id
    INTEGER                                                  :: var_3_id
    INTEGER                                                  :: var_4_id
    INTEGER                                                  :: var_5_id
    INTEGER                                                  :: var_6_id
    INTEGER                                                  :: var_7_id
    INTEGER                                                  :: var_8_id
    INTEGER                                                  :: var_9_id
    INTEGER                                                  :: i, j
    
    ! Creating from the two dimensional lon lat variables, the one dimensional ones:
    DO i = 1, C%NLON
     lon_gcm_1D(i) = lon_gcm(i,1)
    END DO
    DO j = 1, C%NLAT
     lat_gcm_1D(j) = lat_gcm(1,j)
    END DO
    
    ! Write variables to file:
    CALL handle_error(nf90_create(file_gcm, nf90_clobber, ncid))
    ! Define dimensions 
    CALL handle_error(nf90_def_dim(ncid, 'i_lon', C%NLON, lon_dim))
    CALL handle_error(nf90_def_dim(ncid, 'j_lat', C%NLAT, lat_dim))
    ! Define variables
    CALL handle_error(nf90_def_var(ncid, C%gcm_x_axis_name, nf90_float, (/lon_dim/), lon_id_1D))
    CALL handle_error(nf90_def_var(ncid, C%gcm_y_axis_name, nf90_float, (/lat_dim/), lat_id_1D))
    CALL handle_error(nf90_def_var(ncid, field_name_1_gcm , nf90_float, (/lon_dim, lat_dim/), var_1_id))
    IF(PRESENT(field_2_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_2_gcm , nf90_float, (/lon_dim, lat_dim/), var_2_id))
    IF(PRESENT(field_3_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_3_gcm , nf90_float, (/lon_dim, lat_dim/), var_3_id))
    IF(PRESENT(field_4_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_4_gcm , nf90_float, (/lon_dim, lat_dim/), var_4_id))
    IF(PRESENT(field_5_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_5_gcm , nf90_float, (/lon_dim, lat_dim/), var_5_id))
    IF(PRESENT(field_6_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_6_gcm , nf90_float, (/lon_dim, lat_dim/), var_6_id))
    IF(PRESENT(field_7_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_7_gcm , nf90_float, (/lon_dim, lat_dim/), var_7_id))
    IF(PRESENT(field_8_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_8_gcm , nf90_float, (/lon_dim, lat_dim/), var_8_id))
    IF(PRESENT(field_9_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_9_gcm , nf90_float, (/lon_dim, lat_dim/), var_9_id))
    CALL handle_error(nf90_enddef(ncid))
    ! Put variables:
    CALL handle_error(nf90_put_var(ncid, lon_id_1D, lon_gcm_1D, start=(/C%starting_recordnr/)))
    CALL handle_error(nf90_put_var(ncid, lat_id_1D, lat_gcm_1D, start=(/C%starting_recordnr/)))
    CALL handle_error(nf90_put_var(ncid, var_1_id, field_1_gcm,  start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_2_gcm)) CALL handle_error(nf90_put_var(ncid, var_2_id, field_2_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_3_gcm)) CALL handle_error(nf90_put_var(ncid, var_3_id, field_3_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_4_gcm)) CALL handle_error(nf90_put_var(ncid, var_4_id, field_4_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_5_gcm)) CALL handle_error(nf90_put_var(ncid, var_5_id, field_5_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_6_gcm)) CALL handle_error(nf90_put_var(ncid, var_6_id, field_6_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_7_gcm)) CALL handle_error(nf90_put_var(ncid, var_7_id, field_7_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_8_gcm)) CALL handle_error(nf90_put_var(ncid, var_8_id, field_8_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_9_gcm)) CALL handle_error(nf90_put_var(ncid, var_9_id, field_9_gcm, start=(/1, 1, C%starting_recordnr/)))

    ! close the file 
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE write_data_for_gcm



  SUBROUTINE write_data_for_gcm_with_2D_lonlat(file_gcm, lon_gcm, lat_gcm, field_name_1_gcm, field_1_gcm, &
                                                                           field_name_2_gcm, field_2_gcm, &
                                                                           field_name_3_gcm, field_3_gcm, &
                                                                           field_name_4_gcm, field_4_gcm, &
                                                                           field_name_5_gcm, field_5_gcm, &
                                                                           field_name_6_gcm, field_6_gcm, &
                                                                           field_name_7_gcm, field_7_gcm, &
                                                                           field_name_8_gcm, field_8_gcm, &
                                                                           field_name_9_gcm, field_9_gcm)
    USE oblimap_configuration_module, ONLY : dp, C
    USE netcdf
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                   INTENT(IN)           :: file_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lat_gcm
    CHARACTER(LEN=*),                   INTENT(IN)           :: field_name_1_gcm ! name of a 1st GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_2_gcm ! name of a 2nd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_3_gcm ! name of a 3rd GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_4_gcm ! name of a 4th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_5_gcm ! name of a 5th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_6_gcm ! name of a 6th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_7_gcm ! name of a 7th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_8_gcm ! name of a 8th GCM variable
    CHARACTER(LEN=*),                   INTENT(IN), OPTIONAL :: field_name_9_gcm ! name of a 9th GCM variable
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: field_1_gcm      ! values of a 1st GCM field
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_2_gcm      ! values of a 2nd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_3_gcm      ! values of a 3rd GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_4_gcm      ! values of a 4th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_5_gcm      ! values of a 5th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_6_gcm      ! values of a 6th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_7_gcm      ! values of a 7th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_8_gcm      ! values of a 8th GCM field 
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN), OPTIONAL :: field_9_gcm      ! values of a 9th GCM field 
    
    ! Local variables:
    INTEGER                                                  :: ncid
    INTEGER                                                  :: lon_dim
    INTEGER                                                  :: lat_dim
    INTEGER                                                  :: lon_id
    INTEGER                                                  :: lat_id
    INTEGER                                                  :: var_1_id
    INTEGER                                                  :: var_2_id
    INTEGER                                                  :: var_3_id
    INTEGER                                                  :: var_4_id
    INTEGER                                                  :: var_5_id
    INTEGER                                                  :: var_6_id
    INTEGER                                                  :: var_7_id
    INTEGER                                                  :: var_8_id
    INTEGER                                                  :: var_9_id
    
    ! Write variables to file:
    CALL handle_error(nf90_create(file_gcm, nf90_clobber, ncid))
    ! Define dimensions: 
    CALL handle_error(nf90_def_dim(ncid, 'i_lon', C%NLON, lon_dim))
    CALL handle_error(nf90_def_dim(ncid, 'j_lat', C%NLAT, lat_dim))
    ! Define variables:
    CALL handle_error(nf90_def_var(ncid, C%gcm_x_axis_name, nf90_float, (/lon_dim, lat_dim/), lon_id))
    CALL handle_error(nf90_def_var(ncid, C%gcm_y_axis_name, nf90_float, (/lon_dim, lat_dim/), lat_id))
    CALL handle_error(nf90_def_var(ncid, field_name_1_gcm, nf90_float, (/lon_dim, lat_dim/), var_1_id))
    IF(PRESENT(field_2_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_2_gcm , nf90_float, (/lon_dim, lat_dim/), var_2_id))
    IF(PRESENT(field_3_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_3_gcm , nf90_float, (/lon_dim, lat_dim/), var_3_id))
    IF(PRESENT(field_4_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_4_gcm , nf90_float, (/lon_dim, lat_dim/), var_4_id))
    IF(PRESENT(field_5_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_5_gcm , nf90_float, (/lon_dim, lat_dim/), var_5_id))
    IF(PRESENT(field_6_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_6_gcm , nf90_float, (/lon_dim, lat_dim/), var_6_id))
    IF(PRESENT(field_7_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_7_gcm , nf90_float, (/lon_dim, lat_dim/), var_7_id))
    IF(PRESENT(field_8_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_8_gcm , nf90_float, (/lon_dim, lat_dim/), var_8_id))
    IF(PRESENT(field_9_gcm)) CALL handle_error(nf90_def_var(ncid, field_name_9_gcm , nf90_float, (/lon_dim, lat_dim/), var_9_id))
    CALL handle_error(nf90_enddef(ncid))

    ! Put variables:
    CALL handle_error(nf90_put_var(ncid, lon_id, lon_gcm, start=(/1, 1, C%starting_recordnr/)))
    CALL handle_error(nf90_put_var(ncid, lat_id, lat_gcm, start=(/1, 1, C%starting_recordnr/)))
    CALL handle_error(nf90_put_var(ncid, var_1_id, field_1_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_2_gcm)) CALL handle_error(nf90_put_var(ncid, var_2_id, field_2_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_3_gcm)) CALL handle_error(nf90_put_var(ncid, var_3_id, field_3_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_4_gcm)) CALL handle_error(nf90_put_var(ncid, var_4_id, field_4_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_5_gcm)) CALL handle_error(nf90_put_var(ncid, var_5_id, field_5_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_6_gcm)) CALL handle_error(nf90_put_var(ncid, var_6_id, field_6_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_7_gcm)) CALL handle_error(nf90_put_var(ncid, var_7_id, field_7_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_8_gcm)) CALL handle_error(nf90_put_var(ncid, var_8_id, field_8_gcm, start=(/1, 1, C%starting_recordnr/)))
    IF(PRESENT(field_9_gcm)) CALL handle_error(nf90_put_var(ncid, var_9_id, field_9_gcm, start=(/1, 1, C%starting_recordnr/)))

    ! close the file 
    CALL handle_error(nf90_close(ncid))
  END SUBROUTINE write_data_for_gcm_with_2D_lonlat



  SUBROUTINE handle_error(stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE
  
    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF(stat /= nf90_noerr) THEN
     IF(PRESENT(message)) THEN
      WRITE(UNIT=*,FMT='(A,A,A,A)') 'ERROR: netCDF failed (in read_and_write_gcm_im_data_module) because: ', TRIM(nf90_strerror(stat)), ' concerning: ', message
     ELSE
      WRITE(UNIT=*,FMT='(A,A)')     'ERROR: netCDF failed (in read_and_write_gcm_im_data_module) because: ', TRIM(nf90_strerror(stat))
     END IF
     STOP
    END IF
  END SUBROUTINE handle_error
      
END MODULE oblimap_read_and_write_module
