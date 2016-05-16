# oblimap-1.0
oblimap's first release (2010)

Abstract. Here, we present a mapping method OBLIMAP, which projects and interpolates fields like surface temperature, surface mass balance, and surface height between a geographical based coordinate system of a General Circulation Model (GCM) and a rectangular based Ice Model (IM). We derive an oblique stereographic projection and its inverse, which holds for any area at the Earth's surface, and which can be combined with two different interpolation methods. The first one is suited to interpolate the projected fields of a coarse GCM grid on a fine meshed IM grid. The second one is appropriate for the opposite case. Both grids are allowed to be arbitrary and irregularly spaced. Therefore the OBLIMAP technique is suitable for any GCM-IM combination. After a first scan of the GCM grid coordinates and the specification of the IM grid, fast mapping of various fields is possible. To and fro (GCM-IM-GCM) mapping tests with the Climate Community System Model (CCSM) at T42 resolution (~313 km) and the Regional Atmospheric Climate Model (RACMO) at ~11 km and ~55 km, show average temperature differences of less than 0.1 K with small standard deviations. OBLIMAP, available at GMD, is an accurate, robust and well-documented mapping method for coupling an IM with a GCM or to map state of the art initial and forcing fields available at geographical coordinates to any local IM grid with an optimal centered oblique projection. Currently, the oblique stereographic and the oblique Lambert azimuthal equal-area projections for both the sphere and the ellipsoid are implemented in OBLIMAP.

http://www.geosci-model-dev.net/3/13/2010/gmd-3-13-2010.html




To compile the OBLIMAP source:

 Check which fortran compiler you are using. In the Makefile default the
  
  Makefile.pgf90
 
 is included, however also one of below can be included:
  Makefile.g95
  Makefile.gfortran
  Makefile.mpfort
 
 In the Makefile.* one should check the NETCDF path, in our default it is:
  
  NETCDF_DIR = /usr/local
 
 After these checks / addaptions you can compile the OBLIMAP source by
 
  cd src/
  make all
  
 To run:
 
  cd ../
  ./oblimap_to_and_fro_mapping.csh
  
 or with specifying a CONFIG file:
  
  ./oblimap_to_and_fro_mapping.csh CONFIGs/CONFIG_mapped_ccsm_Greenland


Thomas Reerink
Januari 2010
