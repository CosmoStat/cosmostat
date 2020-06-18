
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetCoordsysGrp

FUNCTION PIOGETCOORDSYSGRP,Coordsys,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Coordsys_TMP=BYTARR(128)
    Coordsys_TMP(*)=0
    if (N_ELEMENTS(Coordsys) GT 0) THEN if (STRLEN(Coordsys) GT 0) THEN Coordsys_TMP(0:STRLEN(Coordsys)-1)=BYTE(Coordsys)

 PIOGETCOORDSYSGRP=call_external(PIOLibIDLSO,'piogetcoordsysgrp_tempoidl', $
        Coordsys_TMP, $
        MyGroup, $
               /L64_VALUE) 
  Coordsys=STRING(Coordsys_TMP)

 RETURN,PIOGETCOORDSYSGRP
END
