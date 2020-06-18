
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateVECTGrp

FUNCTION PIOCREATEVECTGRP,Groupname
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)

 PIOCREATEVECTGRP=call_external(PIOLibIDLSO,'piocreatevectgrp_tempoidl', $
        Groupname_TMP, $
               /L64_VALUE) 

 RETURN,PIOCREATEVECTGRP
END
