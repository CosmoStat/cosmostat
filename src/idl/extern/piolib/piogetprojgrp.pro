
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetProjGrp

FUNCTION PIOGETPROJGRP,Proj,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Proj_TMP=BYTARR(128)
    Proj_TMP(*)=0
    if (N_ELEMENTS(Proj) GT 0) THEN if (STRLEN(Proj) GT 0) THEN Proj_TMP(0:STRLEN(Proj)-1)=BYTE(Proj)

 PIOGETPROJGRP=call_external(PIOLibIDLSO,'piogetprojgrp_tempoidl', $
        Proj_TMP, $
        MyGroup, $
               /L64_VALUE) 
  Proj=STRING(Proj_TMP)

 RETURN,PIOGETPROJGRP
END
