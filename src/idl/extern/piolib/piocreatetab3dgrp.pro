
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateTAB3DGrp

FUNCTION PIOCREATETAB3DGRP,Groupname,NAXIS1,NAXIS2
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
    if (N_ELEMENTS(NAXIS1) EQ 0) THEN     NAXIS1=LONG64(0) ELSE    NAXIS1=LONG64(NAXIS1)
     if (N_ELEMENTS(NAXIS2) EQ 0) THEN     NAXIS2=LONG64(0) ELSE    NAXIS2=LONG64(NAXIS2)
 
 PIOCREATETAB3DGRP=call_external(PIOLibIDLSO,'piocreatetab3dgrp_tempoidl', $
        Groupname_TMP, $
        NAXIS1, $
        NAXIS2, $
               /L64_VALUE) 

 RETURN,PIOCREATETAB3DGRP
END
