
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateTAB2DGrp

FUNCTION PIOCREATETAB2DGRP,Groupname,NumberofColumn
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
    if (N_ELEMENTS(NumberofColumn) EQ 0) THEN     NumberofColumn=LONG64(0) ELSE    NumberofColumn=LONG64(NumberofColumn)
 
 PIOCREATETAB2DGRP=call_external(PIOLibIDLSO,'piocreatetab2dgrp_tempoidl', $
        Groupname_TMP, $
        NumberofColumn, $
               /L64_VALUE) 

 RETURN,PIOCREATETAB2DGRP
END
