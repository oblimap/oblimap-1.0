! File name: oblimap_mapping_module.f90          
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

MODULE oblimap_mapping_module

CONTAINS
  SUBROUTINE fast_mapping(fast_input_file_name, N_row_input, N_column_input, N_row_mapped, N_column_mapped, input_field_1, mapped_field_1, &
                                                                                                            input_field_2, mapped_field_2, &
                                                                                                            input_field_3, mapped_field_3, &
                                                                                                            input_field_4, mapped_field_4, &
                                                                                                            input_field_5, mapped_field_5, &
                                                                                                            input_field_6, mapped_field_6, &
                                                                                                            input_field_7, mapped_field_7, &
                                                                                                            input_field_8, mapped_field_8, &
                                                                                                            input_field_9, mapped_field_9)
    ! With this routine the fast mapping can be done. 
    !   mapping = (inverse) oblique projection + interpolation
    ! This routine reads the 'fast mapping' input file. It can handle the four situations of 
    !  oblique projection + quadrant interpolation
    !  oblique projection + radius   interpolation
    !  inverse oblique projection + quadrant interpolation
    !  inverse oblique projection + radius   interpolation
    USE oblimap_configuration_module, ONLY : dp, C
    IMPLICIT NONE

    ! Input variable:
    CHARACTER(LEN=*),                                   INTENT(IN)            :: fast_input_file_name 
    INTEGER,                                            INTENT(IN)            :: N_row_input
    INTEGER,                                            INTENT(IN)            :: N_column_input
    INTEGER,                                            INTENT(IN)            :: N_row_mapped
    INTEGER,                                            INTENT(IN)            :: N_column_mapped
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN)            :: input_field_1
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_2
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_3
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_4
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_5
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_6
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_7
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_8
    REAL(dp), DIMENSION(N_row_input,  N_column_input),  INTENT(IN),  OPTIONAL :: input_field_9

    ! Output variables:
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT)           :: mapped_field_1
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_2
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_3
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_4
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_5
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_6
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_7
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_8
    REAL(dp), DIMENSION(N_row_mapped, N_column_mapped), INTENT(OUT), OPTIONAL :: mapped_field_9
           
    ! Local variables:
    INTEGER                                                                   :: passing_header_line  ! Counter which walks over the header lines
    INTEGER                                                                   :: total_mapped_points  ! Amount of affected/mapped/target grid points by the mapping
    INTEGER                                                                   :: count_mapped_points  ! Counter which counts over the affected/mapped/target points
    INTEGER                                                                   :: row_mapped           ! row index of the affected/mapped/target points
    INTEGER                                                                   :: column_mapped        ! column index of the affected/mapped/target points
    INTEGER                                                                   :: total_contributions  ! Number contributing points used to estimate the field value for the grid point
    INTEGER                                                                   :: contribution         ! Counter for the involved IM points for each GCM point
    INTEGER                                                                   :: row_input            ! row index of the input points
    INTEGER                                                                   :: column_input         ! column index of the input points
    REAL(dp)                                                                  :: distance             ! distance over the projection plane between projected point and grid point
    REAL(dp)                                                                  :: numerator_1          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_2          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_3          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_4          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_5          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_6          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_7          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_8          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: numerator_9          ! The   numerator in equation (2.19) in Reerink et al. (2010), the Shepard formula.
    REAL(dp)                                                                  :: denumerator          ! The denumerator in equation (2.19) in Reerink et al. (2010), the Shepard formula. 
    CHARACTER(1)                                                              :: end_of_line         

    ! Opening the 'fast input' file:
    OPEN(UNIT=117, FILE=TRIM(fast_input_file_name))

    ! Ignoring the header while reading the header:
    DO passing_header_line = 1, 11
     READ(UNIT=117, FMT='(A)') end_of_line
    END DO

    ! Reading the total number of mapped points. Each line in the file contains the necessary information of all the 
    ! contributing points for one target grid point. The number of mapped (or target) points equals the number lines 
    READ(UNIT=117, FMT='(I10)') total_mapped_points

    ! See equation (2.17) and equation (2.19) in Reerink et al. (2010), both cases are treated with the same code:  
    DO count_mapped_points = 1, total_mapped_points
      numerator_1  = 0._dp
      IF(PRESENT(mapped_field_2)) numerator_2  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_3  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_4  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_5  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_6  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_7  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_8  = 0._dp
      IF(PRESENT(mapped_field_3)) numerator_9  = 0._dp
      denumerator  = 0._dp
      READ(UNIT=117, FMT='(3I6)', ADVANCE='NO') row_mapped, column_mapped, total_contributions
      DO contribution = 1, total_contributions
        READ(UNIT=117, FMT='(2I6,E23.15)', ADVANCE='NO') row_input, column_input, distance

        numerator_1 = numerator_1 + input_field_1(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_2)) numerator_2 = numerator_2 + input_field_2(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_3)) numerator_3 = numerator_3 + input_field_3(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_4)) numerator_4 = numerator_4 + input_field_4(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_5)) numerator_5 = numerator_5 + input_field_5(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_6)) numerator_6 = numerator_6 + input_field_6(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_7)) numerator_7 = numerator_7 + input_field_7(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_8)) numerator_8 = numerator_8 + input_field_8(row_input,column_input) / (distance**C%shephard_exponent_radius)
        IF(PRESENT(mapped_field_9)) numerator_9 = numerator_9 + input_field_9(row_input,column_input) / (distance**C%shephard_exponent_radius)
        denumerator  = denumerator + (1._dp / (distance**C%shephard_exponent_radius))
      END DO
      READ(UNIT=117, FMT='(A)') end_of_line
      IF(total_contributions > 0) THEN
       mapped_field_1(row_mapped,column_mapped) = numerator_1 / denumerator
       IF(PRESENT(mapped_field_2)) mapped_field_2(row_mapped,column_mapped) = numerator_2 / denumerator
       IF(PRESENT(mapped_field_3)) mapped_field_3(row_mapped,column_mapped) = numerator_3 / denumerator
       IF(PRESENT(mapped_field_4)) mapped_field_4(row_mapped,column_mapped) = numerator_4 / denumerator
       IF(PRESENT(mapped_field_5)) mapped_field_5(row_mapped,column_mapped) = numerator_5 / denumerator
       IF(PRESENT(mapped_field_6)) mapped_field_6(row_mapped,column_mapped) = numerator_6 / denumerator
       IF(PRESENT(mapped_field_7)) mapped_field_7(row_mapped,column_mapped) = numerator_7 / denumerator
       IF(PRESENT(mapped_field_8)) mapped_field_8(row_mapped,column_mapped) = numerator_8 / denumerator
       IF(PRESENT(mapped_field_9)) mapped_field_9(row_mapped,column_mapped) = numerator_9 / denumerator
      END IF
    END DO
    
    ! Closing the 'fast input' file:
    CLOSE(UNIT=117)
  END SUBROUTINE fast_mapping
      
END MODULE oblimap_mapping_module
