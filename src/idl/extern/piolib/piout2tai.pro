
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOUT2TAI

FUNCTION PIOUT2TAI,UT,TAI
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    UT_TMP=BYTARR(128)
    UT_TMP(*)=0
    if (N_ELEMENTS(UT) GT 0) THEN if (STRLEN(UT) GT 0) THEN UT_TMP(0:STRLEN(UT)-1)=BYTE(UT)
    if (N_ELEMENTS(TAI) EQ 0) THEN     TAI=DOUBLE(0) ELSE    TAI=DOUBLE(TAI)
 
 PIOUT2TAI=call_external(PIOLibIDLSO,'piout2tai_tempoidl', $
        UT_TMP, $
        TAI, $
               /L64_VALUE) 

 RETURN,PIOUT2TAI
END
