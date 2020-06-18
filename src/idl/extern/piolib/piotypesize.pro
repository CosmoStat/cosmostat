
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOTypeSize

FUNCTION PIOTYPESIZE,typeIn
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    typeIn_TMP=BYTARR(128)
    typeIn_TMP(*)=0
    if (N_ELEMENTS(typeIn) GT 0) THEN if (STRLEN(typeIn) GT 0) THEN typeIn_TMP(0:STRLEN(typeIn)-1)=BYTE(typeIn)

 PIOTYPESIZE=call_external(PIOLibIDLSO,'piotypesize_tempoidl', $
        typeIn_TMP, $
               /L64_VALUE) 

 RETURN,PIOTYPESIZE
END
