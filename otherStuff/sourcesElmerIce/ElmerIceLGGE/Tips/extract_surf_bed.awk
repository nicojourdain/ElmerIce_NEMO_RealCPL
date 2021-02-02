awk '{if ($3==1) printf "%+8.8e %+8.8e\n", $24, $25}' surfMISMIP3a_s01.dat > mismip_bed.dat 
awk '{if ($3==3) printf "%+8.8e %+8.8e\n", $24, $25}' surfMISMIP3a_s01.dat > mismip_surf.dat 
