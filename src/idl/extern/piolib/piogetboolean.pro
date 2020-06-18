
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetBoolean

FUNCTION PIOGETBOOLEAN,ParamName,MyPar
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ParamName_TMP=BYTARR(128)
    ParamName_TMP(*)=0
    if (N_ELEMENTS(ParamName) GT 0) THEN if (STRLEN(ParamName) GT 0) THEN ParamName_TMP(0:STRLEN(ParamName)-1)=BYTE(ParamName)
 IF (N_ELEMENTS(MyPar) EQ 0) THEN MyPar=LONG64(0) ELSE MyPar=LONG64(MyPar)

 PIOGETBOOLEAN=call_external(PIOLibIDLSO,'piogetboolean_tempoidl', $
        ParamName_TMP, $
        MyPar, $
               /L64_VALUE) 

 RETURN,PIOGETBOOLEAN
END
