
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOpenParFile

FUNCTION PIOOPENPARFILE,ParName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ParName_TMP=BYTARR(128)
    ParName_TMP(*)=0
    if (N_ELEMENTS(ParName) GT 0) THEN if (STRLEN(ParName) GT 0) THEN ParName_TMP(0:STRLEN(ParName)-1)=BYTE(ParName)

 PIOOPENPARFILE=call_external(PIOLibIDLSO,'pioopenparfile_tempoidl', $
        ParName_TMP, $
               /L64_VALUE) 

 RETURN,PIOOPENPARFILE
END
