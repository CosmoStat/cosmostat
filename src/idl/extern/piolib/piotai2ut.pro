
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOTAI2UT

FUNCTION PIOTAI2UT,TAI,UT
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(TAI) EQ 0) THEN     TAI=DOUBLE(0) ELSE    TAI=DOUBLE(TAI)
     UT_TMP=BYTARR(128)
    UT_TMP(*)=0
    if (N_ELEMENTS(UT) GT 0) THEN if (STRLEN(UT) GT 0) THEN UT_TMP(0:STRLEN(UT)-1)=BYTE(UT)

 PIOTAI2UT=call_external(PIOLibIDLSO,'piotai2ut_tempoidl', $
        TAI, $
        UT_TMP, $
               /L64_VALUE) 

 RETURN,PIOTAI2UT
END
