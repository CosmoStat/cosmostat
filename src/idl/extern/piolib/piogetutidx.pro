
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetUTIdx

FUNCTION PIOGETUTIDX,GrpName,UTC
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GrpName_TMP=BYTARR(128)
    GrpName_TMP(*)=0
    if (N_ELEMENTS(GrpName) GT 0) THEN if (STRLEN(GrpName) GT 0) THEN GrpName_TMP(0:STRLEN(GrpName)-1)=BYTE(GrpName)
    UTC_TMP=BYTARR(128)
    UTC_TMP(*)=0
    if (N_ELEMENTS(UTC) GT 0) THEN if (STRLEN(UTC) GT 0) THEN UTC_TMP(0:STRLEN(UTC)-1)=BYTE(UTC)

 PIOGETUTIDX=call_external(PIOLibIDLSO,'piogetutidx_tempoidl', $
        GrpName_TMP, $
        UTC_TMP, $
               /L64_VALUE) 

 RETURN,PIOGETUTIDX
END
