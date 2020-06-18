
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateTIMESTAMPGrp

FUNCTION PIOCREATETIMESTAMPGRP,Groupname,list_TRII
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
    list_TRII_TMP=BYTARR(128)
    list_TRII_TMP(*)=0
    if (N_ELEMENTS(list_TRII) GT 0) THEN if (STRLEN(list_TRII) GT 0) THEN list_TRII_TMP(0:STRLEN(list_TRII)-1)=BYTE(list_TRII)

 PIOCREATETIMESTAMPGRP=call_external(PIOLibIDLSO,'piocreatetimestampgrp_tempoidl', $
        Groupname_TMP, $
        list_TRII_TMP, $
               /L64_VALUE) 

 RETURN,PIOCREATETIMESTAMPGRP
END
