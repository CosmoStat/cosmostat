
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetListSz

FUNCTION PIOGETLISTSZ,ListName,MyPar
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ListName_TMP=BYTARR(128)
    ListName_TMP(*)=0
    if (N_ELEMENTS(ListName) GT 0) THEN if (STRLEN(ListName) GT 0) THEN ListName_TMP(0:STRLEN(ListName)-1)=BYTE(ListName)
 IF (N_ELEMENTS(MyPar) EQ 0) THEN MyPar=LONG64(0) ELSE MyPar=LONG64(MyPar)

 PIOGETLISTSZ=call_external(PIOLibIDLSO,'piogetlistsz_tempoidl', $
        ListName_TMP, $
        MyPar, $
               /L64_VALUE) 

 RETURN,PIOGETLISTSZ
END
