
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateROIGrp

FUNCTION PIOCREATEROIGRP,Groupname,Init_Size
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
IF (N_PARAMS() LT 2) THEN BEGIN
    Init_Size=long64(0)
END ELSE BEGIN
    Init_Size=long64(Init_Size)
END

 PIOCREATEROIGRP=call_external(PIOLibIDLSO,'piocreateroigrp_tempoidl', $
        Groupname_TMP, $
        Init_Size, $
               /L64_VALUE) 

 RETURN,PIOCREATEROIGRP
END
