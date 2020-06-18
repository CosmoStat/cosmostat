
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTAB2DLines

FUNCTION PIOGETTAB2DLINES,FirstLine,LastLine,ObjectName,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(FirstLine) EQ 0) THEN     FirstLine=LONG64(0) ELSE    FirstLine=LONG64(FirstLine)
     if (N_ELEMENTS(LastLine) EQ 0) THEN     LastLine=LONG64(0) ELSE    LastLine=LONG64(LastLine)
     ObjectName_TMP=BYTARR(128)
    ObjectName_TMP(*)=0
    if (N_ELEMENTS(ObjectName) GT 0) THEN if (STRLEN(ObjectName) GT 0) THEN ObjectName_TMP(0:STRLEN(ObjectName)-1)=BYTE(ObjectName)

 PIOGETTAB2DLINES=call_external(PIOLibIDLSO,'piogettab2dlines_tempoidl', $
        FirstLine, $
        LastLine, $
        ObjectName_TMP, $
        MyGroup, $
               /L64_VALUE) 
  ObjectName=STRING(ObjectName_TMP)

 RETURN,PIOGETTAB2DLINES
END
