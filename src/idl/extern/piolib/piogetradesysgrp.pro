
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetRadesysGrp

FUNCTION PIOGETRADESYSGRP,RadSysN,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    RadSysN_TMP=BYTARR(128)
    RadSysN_TMP(*)=0
    if (N_ELEMENTS(RadSysN) GT 0) THEN if (STRLEN(RadSysN) GT 0) THEN RadSysN_TMP(0:STRLEN(RadSysN)-1)=BYTE(RadSysN)

 PIOGETRADESYSGRP=call_external(PIOLibIDLSO,'piogetradesysgrp_tempoidl', $
        RadSysN_TMP, $
        MyGroup, $
               /L64_VALUE) 
  RadSysN=STRING(RadSysN_TMP)

 RETURN,PIOGETRADESYSGRP
END
