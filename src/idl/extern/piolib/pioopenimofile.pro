
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOpenIMOFile

FUNCTION PIOOPENIMOFILE,File,openm
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    File_TMP=BYTARR(128)
    File_TMP(*)=0
    if (N_ELEMENTS(File) GT 0) THEN if (STRLEN(File) GT 0) THEN File_TMP(0:STRLEN(File)-1)=BYTE(File)
    openm_TMP=BYTARR(128)
    openm_TMP(*)=0
    if (N_ELEMENTS(openm) GT 0) THEN if (STRLEN(openm) GT 0) THEN openm_TMP(0:STRLEN(openm)-1)=BYTE(openm)

 PIOOPENIMOFILE=call_external(PIOLibIDLSO,'pioopenimofile_tempoidl', $
        File_TMP, $
        openm_TMP, $
               /L64_VALUE) 

 RETURN,PIOOPENIMOFILE
END
