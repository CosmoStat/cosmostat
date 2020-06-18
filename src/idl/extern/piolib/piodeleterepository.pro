
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteRepository

FUNCTION PIODELETEREPOSITORY,RepositoryName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    RepositoryName_TMP=BYTARR(128)
    RepositoryName_TMP(*)=0
    if (N_ELEMENTS(RepositoryName) GT 0) THEN if (STRLEN(RepositoryName) GT 0) THEN RepositoryName_TMP(0:STRLEN(RepositoryName)-1)=BYTE(RepositoryName)

 PIODELETEREPOSITORY=call_external(PIOLibIDLSO,'piodeleterepository_tempoidl', $
        RepositoryName_TMP, $
               /L64_VALUE) 

 RETURN,PIODELETEREPOSITORY
END
