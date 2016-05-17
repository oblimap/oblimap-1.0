! File name: oblimap_projection_module.f90          
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

MODULE oblimap_projection_module
    
CONTAINS
  SUBROUTINE oblique_sg_projection(lambda, phi, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate 
    ! system, with coordinates (x,y).
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: t_P_prime
    
    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%deg2rad * phi       
    lambda_P = C%deg2rad * lambda

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(C%alpha_stereographic)) / (1._dp + COS(phi_P) * COS(C%phi_M) * COS(lambda_P - C%lambda_M) + SIN(phi_P) * SIN(C%phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  C%R_earth * (COS(phi_P) * SIN(lambda_P - C%lambda_M)) * t_P_prime
    y_IM_P_prime =  C%R_earth * (SIN(phi_P) * COS(C%phi_M) - (COS(phi_P) * SIN(C%phi_M)) * COS(lambda_P - C%lambda_M)) * t_P_prime
  END SUBROUTINE oblique_sg_projection



  SUBROUTINE inverse_oblique_sg_projection(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the 
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime
    REAL(dp)              :: y_3D_P_prime
    REAL(dp)              :: z_3D_P_prime
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P
    REAL(dp)              :: y_3D_P
    REAL(dp)              :: z_3D_P
    
    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = C%R_earth * COS(C%alpha_stereographic) * COS(C%lambda_M) * COS(C%phi_M) - SIN(C%lambda_M) * x_IM_P_prime - COS(C%lambda_M) * SIN(C%phi_M) * y_IM_P_prime
    y_3D_P_prime = C%R_earth * COS(C%alpha_stereographic) * SIN(C%lambda_M) * COS(C%phi_M) + COS(C%lambda_M) * x_IM_P_prime - SIN(C%lambda_M) * SIN(C%phi_M) * y_IM_P_prime
    z_3D_P_prime = C%R_earth * COS(C%alpha_stereographic) *                   SIN(C%phi_M)                                  +                   COS(C%phi_M) * y_IM_P_prime
    
    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(C%lambda_M) * COS(C%phi_M) * x_3D_P_prime  +  SIN(C%lambda_M) * COS(C%phi_M) * y_3D_P_prime  +  SIN(C%phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * C%R_earth**2 + 2._dp * C%R_earth * a) / (C%R_earth**2 + 2._dp * C%R_earth * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  C%R_earth * COS(C%lambda_M) * COS(C%phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  C%R_earth * SIN(C%lambda_M) * COS(C%phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  C%R_earth *                   SIN(C%phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
     lambda_P = 180._dp + C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
     lambda_P =           C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 360._dp + C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
     lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
     lambda_P =   0._dp
    END IF

   ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
   IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
    phi_P = C%rad2deg * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2)) 
   ELSE IF(z_3D_P >  0._dp) THEN
    phi_P =   90._dp
   ELSE IF(z_3D_P <  0._dp) THEN
    phi_P =  -90._dp
   END IF
  END SUBROUTINE inverse_oblique_sg_projection



  SUBROUTINE inverse_oblique_sg_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection the 
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C     ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (20-18) on page 159 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
    ! See equation (21-15) on page 159 Snyder (1987), because the denumerator is always positive this ATAN doesn't 
    ! need a correction like note 2 on page ix in Snyder (1987):
    angle_C  = 2._dp * ATAN(rho / ((1._dp + COS(C%alpha_stereographic)) * C%R_earth))
    
    ! See equation (20-14) on page 158 Snyder (1987):
    phi_P    = C%rad2deg * ( ASIN(COS(angle_C) * SIN(C%phi_M) + ((y_IM_P_prime * SIN(angle_C) * COS(C%phi_M)) / rho)) )
    
    ! See equation (20-15) on page 159 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = rho * COS(C%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%phi_M) * SIN(angle_C)
    lambda_P    = C%rad2deg * (C%lambda_M + arctanges_quotient(numerator, denumerator))
    
    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    
    ! In case point P coincides with M (see condition at the first line of page  159 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%rad2deg * C%lambda_M
     phi_P    = C%rad2deg * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_sg_projection_snyder



  SUBROUTINE oblique_laea_projection_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection the 
    ! longitude-latitude coordinates which coincide with the GCM grid points to the rectangular IM 
    ! coordinate system, with coordinates (x,y).
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: t_P_prime
    
    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%deg2rad * phi       
    lambda_P = C%deg2rad * lambda

    ! See equation (21-4) on page 185 of Snyder (1987):
    t_P_prime = SQRT(2._dp / (1._dp + COS(phi_P) * COS(C%phi_M) * COS(lambda_P - C%lambda_M) + SIN(phi_P) * SIN(C%phi_M)))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010), page 185 of Snyder (1987):
    x_IM_P_prime =  C%R_earth * (COS(phi_P) * SIN(lambda_P - C%lambda_M)) * t_P_prime
    y_IM_P_prime =  C%R_earth * (SIN(phi_P) * COS(C%phi_M) - (COS(phi_P) * SIN(C%phi_M)) * COS(lambda_P - C%lambda_M)) * t_P_prime
  END SUBROUTINE oblique_laea_projection_snyder



  SUBROUTINE inverse_oblique_laea_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection 
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C              ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (20-18) on page 187 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
    ! See equation (24-16) on page 187 Snyder (1987):
    angle_C  = 2._dp * ASIN(rho / (2._dp * C%R_earth))
    
    ! See equation (20-14) on page 186 Snyder (1987):
    phi_P    = C%rad2deg * ( ASIN(COS(angle_C) * SIN(C%phi_M) + ((y_IM_P_prime * SIN(angle_C) * COS(C%phi_M)) / rho)) )
    
    ! See equation (20-15) on page 186 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = rho * COS(C%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%phi_M) * SIN(angle_C)
    lambda_P    = C%rad2deg * (C%lambda_M + arctanges_quotient(numerator, denumerator))
    
    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    
    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%rad2deg * C%lambda_M
     phi_P    = C%rad2deg * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_laea_projection_snyder



  SUBROUTINE oblique_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with Snyder's oblique stereographic projection for the ellipsoid
    ! the the longitude-latitude coordinates which coincide with the GCM grid points to 
    ! the rectangular IM coordinate system, with coordinates (x,y).
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P    ! phi in Snyder
    REAL(dp)              :: lambda_P ! lambda in Snyder
    REAL(dp)              :: chi_P    ! chi in Snyder
    REAL(dp)              :: A
    
    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%deg2rad * phi       
    lambda_P = C%deg2rad * lambda

    ! See equations (3-1a) and (21-27) on page 160 in Snyder (1987):
    chi_P = 2._dp * ATAN(SQRT(((1._dp +       SIN(phi_P)) / (1._dp -       SIN(phi_P))) * &
                              ((1._dp - C%e * SIN(phi_P)) / (1._dp + C%e * SIN(phi_P)))**(C%e))) - 0.5_dp * C%pi
    A     = C%akm / (COS(C%chi_M) * (1._dp + SIN(C%chi_M) * SIN(chi_P) + COS(C%chi_M) * COS(chi_P) * COS(lambda_P - C%lambda_M)))


    ! See equations (21-24) and (21-25) on page 160 in Snyder (1987):
    x_IM_P_prime =  A * COS(chi_P) * SIN(lambda_P - C%lambda_M)
    y_IM_P_prime =  A * (COS(C%chi_M) * SIN(chi_P) - SIN(C%chi_M) * COS(chi_P) * COS(lambda_P - C%lambda_M))
  END SUBROUTINE oblique_sg_projection_ellipsoid_snyder



  SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection for the ellipsoid 
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C     ! In radians
    REAL(dp)              :: chi_P       ! chi in Snyder
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (20-18) on page 162 Snyder (1987):
    rho     = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
    ! See equation (21-38) on page 162 Snyder (1987):
    angle_C = 2._dp * ATAN(rho * COS(C%chi_M) / C%akm)
    
    ! See equations (21-37) on page 161 in Snyder (1987):
    chi_P   = ASIN(COS(angle_C) * SIN(C%chi_M) + y_IM_P_prime * SIN(angle_C) * COS(C%chi_M) / rho)

    ! See equation (3-5) on page 162 instead of equation (3-4) on page 161 Snyder (1987):
    phi_P = C%rad2deg * (chi_P + &
            (C%e**2 / 2._dp + 5._dp * C%e**4 / 24._dp +          C%e**6 / 12._dp  +   13._dp * C%e**8 /    360._dp) * SIN(2._dp * chi_P) + &
            (                 7._dp * C%e**4 / 48._dp + 29._dp * C%e**6 / 240._dp +  811._dp * C%e**8 /  11520._dp) * SIN(4._dp * chi_P) + &
            (                                            7._dp * C%e**6 / 120._dp +   81._dp * C%e**8 /   1120._dp) * SIN(6._dp * chi_P) + &
            (                                                                       4279._dp * C%e**8 / 161280._dp) * SIN(8._dp * chi_P)) 
    
    ! See equation (21-36) on page 161 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = rho * COS(C%chi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%chi_M) * SIN(angle_C)
    lambda_P    = C%rad2deg * (C%lambda_M + arctanges_quotient(numerator, denumerator))
    
    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    
    ! In case point P coincides with M (see condition at the first line of page  159 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%rad2deg * C%lambda_M
     phi_P    = C%rad2deg * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder



  SUBROUTINE oblique_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection for 
    ! the ellipsoid the longitude-latitude coordinates which coincide with the GCM grid points to 
    ! the rectangular IM coordinate system, with coordinates (x,y).
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: q_P      ! q in Snyder
    REAL(dp)              :: beta
    REAL(dp)              :: B
    
    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%deg2rad * phi       
    lambda_P = C%deg2rad * lambda
    
    ! See equation (3-12) on page 187 in Snyder (1987):
    q_P = (1._dp - C%e**2) * ((SIN(phi_P) / (1._dp - (C%e * SIN(phi_P))**2)) - (1._dp / (2._dp * C%e)) * LOG((1._dp - C%e * SIN(phi_P)) / (1._dp + C%e * SIN(phi_P)))) 
    ! See equation (3-11) on page 187 in Snyder (1987):
    beta = ASIN(q_P / C%q_polar)
    ! See equation (24-19) on page 187 in Snyder (1987):
    B = C%R_q_polar * SQRT(2._dp / (1._dp + SIN(C%beta_M) * SIN(beta) + COS(C%beta_M) * COS(beta) * COS(lambda_P - C%lambda_M)))

    ! See equation (24-17) and (24-18) on page 187 in Snyder (1987):
    x_IM_P_prime = B * C%D * COS(beta) * SIN(lambda_P - C%lambda_M)
    y_IM_P_prime = (B / C%D) * (COS(C%beta_M) * SIN(beta) - SIN(C%beta_M) * COS(beta) * COS(lambda_P - C%lambda_M))
  END SUBROUTINE oblique_laea_projection_ellipsoid_snyder



  SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection for 
    ! the ellipsoid the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C           ! In radians
    REAL(dp)              :: beta              ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (24-28) on page 189 Snyder (1987):
    rho      = SQRT((x_IM_P_prime / C%D)**2 + (y_IM_P_prime / C%D)**2)
    ! See equation (24-29) on page 189 Snyder (1987):
    angle_C  = 2._dp * ASIN(rho / (2._dp * C%R_q_polar))
    
    ! See equation (24-30) on page 189 Snyder (1987):
    beta = ASIN(COS(angle_C) * SIN(C%beta_M) + (C%D * y_IM_P_prime * SIN(angle_C) * COS(C%beta_M) / rho))
    
    ! See equation (3-18) on page 189 instead of equation (3-16) on page 188 Snyder (1987):
    phi_P = C%rad2deg * (beta + &
            (C%e**2 / 3._dp + 31._dp * C%e**4 / 180._dp + 517._dp * C%e**6 /  5040._dp) * SIN(2._dp * beta) + &
            (                 23._dp * C%e**4 / 360._dp + 251._dp * C%e**6 /  3780._dp) * SIN(4._dp * beta) + &
            (                                             761._dp * C%e**6 / 45360._dp) * SIN(6._dp * beta)) 
    
    ! See equation (20-26) on page 188 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = C%D * rho * COS(C%beta_M) * COS(angle_C) - C%D**2 * y_IM_P_prime * SIN(C%beta_M) * SIN(angle_C)
    lambda_P    = C%rad2deg * (C%lambda_M + arctanges_quotient(numerator, denumerator))
    
    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    
    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%rad2deg * C%lambda_M
     phi_P    = C%rad2deg * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder



  FUNCTION arctanges_quotient(numerator, denumerator) RESULT(angle)
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: numerator  
    REAL(dp), INTENT(IN)  :: denumerator

    ! Result variable:
    REAL(dp)              :: angle       ! In radians

    ! Local variables:
    REAL(dp)              :: quadrant_correction
    
    ! See note 2 on page ix in Snyder (1987), to distinguish between the quadrants:
    quadrant_correction = 0._dp
    IF(denumerator <  0._dp) quadrant_correction =          C%pi
    IF(denumerator == 0._dp) quadrant_correction = 0.5_dp * C%pi
    IF(numerator   <  0._dp) quadrant_correction = - quadrant_correction
    
    angle = ATAN(numerator / denumerator) + quadrant_correction
  END FUNCTION arctanges_quotient

END MODULE oblimap_projection_module
