
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCoordsysGrp

FUNCTION PIOCOORDSYSGRP,coordsys,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    coordsys_TMP=BYTARR(128)
    coordsys_TMP(*)=0
    if (N_ELEMENTS(coordsys) GT 0) THEN if (STRLEN(coordsys) GT 0) THEN coordsys_TMP(0:STRLEN(coordsys)-1)=BYTE(coordsys)

 PIOCOORDSYSGRP=call_external(PIOLibIDLSO,'piocoordsysgrp_tempoidl', $
        coordsys_TMP, $
        MyGroup, $
               /L64_VALUE) 
  coordsys=STRING(coordsys_TMP)

 RETURN,PIOCOORDSYSGRP
END
