
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOpenTIMESTAMPGrp

FUNCTION PIOOPENTIMESTAMPGRP,GroupName,mode
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)
    mode_TMP=BYTARR(128)
    mode_TMP(*)=0
    if (N_ELEMENTS(mode) GT 0) THEN if (STRLEN(mode) GT 0) THEN mode_TMP(0:STRLEN(mode)-1)=BYTE(mode)

 PIOOPENTIMESTAMPGRP=call_external(PIOLibIDLSO,'pioopentimestampgrp_tempoidl', $
        GroupName_TMP, $
        mode_TMP, $
               /L64_VALUE) 

 RETURN,PIOOPENTIMESTAMPGRP
END
