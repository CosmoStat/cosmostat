
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOOpenPOIObject

FUNCTION PIOOPENPOIOBJECT,MyObjectName,typein,Command,step,mode,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyObjectName_TMP=BYTARR(128)
    MyObjectName_TMP(*)=0
    if (N_ELEMENTS(MyObjectName) GT 0) THEN if (STRLEN(MyObjectName) GT 0) THEN MyObjectName_TMP(0:STRLEN(MyObjectName)-1)=BYTE(MyObjectName)
    typein_TMP=BYTARR(128)
    typein_TMP(*)=0
    if (N_ELEMENTS(typein) GT 0) THEN if (STRLEN(typein) GT 0) THEN typein_TMP(0:STRLEN(typein)-1)=BYTE(typein)
    Command_TMP=BYTARR(128)
    Command_TMP(*)=0
    if (N_ELEMENTS(Command) GT 0) THEN if (STRLEN(Command) GT 0) THEN Command_TMP(0:STRLEN(Command)-1)=BYTE(Command)
    if (N_ELEMENTS(step) EQ 0) THEN     step=LONG64(0) ELSE    step=LONG64(step)
     mode_TMP=BYTARR(128)
    mode_TMP(*)=0
    if (N_ELEMENTS(mode) GT 0) THEN if (STRLEN(mode) GT 0) THEN mode_TMP(0:STRLEN(mode)-1)=BYTE(mode)
IF (N_PARAMS() LT 6) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOOPENPOIOBJECT=call_external(PIOLibIDLSO,'pioopenpoiobject_tempoidl', $
        MyObjectName_TMP, $
        typein_TMP, $
        Command_TMP, $
        step, $
        mode_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOOPENPOIOBJECT
END
