#! /bin/csh -f
# Thomas Reerink

if($#argv == 0 || $#argv == 1) then
 
 if($#argv == 1) then
  set CONFIG_mapped = $1
 else
  set CONFIG_mapped = CONFIGs/CONFIG_mapped_ccsm_Greenland
 endif

 mkdir -p Input_fast_mapping

 ./src/oblimap_gcm_to_im_program $CONFIG_mapped
 ./src/oblimap_im_to_gcm_program $CONFIG_mapped

else
 echo ' One argument optional, e.g.:'
 echo '  ./src/oblimap_gcm_to_im_program CONFIGs/CONFIG_mapped_ccsm_Greenland'
endif
