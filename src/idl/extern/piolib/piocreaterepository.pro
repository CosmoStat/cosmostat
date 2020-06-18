
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreateRepository

FUNCTION PIOCREATEREPOSITORY,DataRepository,HTMLRepository
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    DataRepository_TMP=BYTARR(128)
    DataRepository_TMP(*)=0
    if (N_ELEMENTS(DataRepository) GT 0) THEN if (STRLEN(DataRepository) GT 0) THEN DataRepository_TMP(0:STRLEN(DataRepository)-1)=BYTE(DataRepository)
    HTMLRepository_TMP=BYTARR(128)
    HTMLRepository_TMP(*)=0
    if (N_ELEMENTS(HTMLRepository) GT 0) THEN if (STRLEN(HTMLRepository) GT 0) THEN HTMLRepository_TMP(0:STRLEN(HTMLRepository)-1)=BYTE(HTMLRepository)

 PIOCREATEREPOSITORY=call_external(PIOLibIDLSO,'piocreaterepository_tempoidl', $
        DataRepository_TMP, $
        HTMLRepository_TMP, $
               /L64_VALUE) 

 RETURN,PIOCREATEREPOSITORY
END
