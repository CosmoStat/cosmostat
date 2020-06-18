
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOWriteVECTObjectMask

FUNCTION PIOWRITEVECTOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 CASE (STRING(TYPE)) OF
"PIOLONG": data = LONG64(data)
"PIOFLAG": data = BYTE(data)
"PIOBYTE": data = BYTE(data)
"PIOINT": data = LONG(data)
"PIOSTRING":BEGIN
    data=PIOSTRING2C(data)
     NBDATA=N_ELEMENTS(data)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(data(I))-1)=BYTE(data(I))
   data=RES
END 
"PIOSHORT": data = FIX(data)
 "PIOFLOAT": data = FLOAT(data)
 "PIODOUBLE": data = DOUBLE(data)
 "PIOCOMPLEX": data = COMPLEX(data)
 "PIODBLCOMPLEX": data = DBLCOMPLEX(data)
ELSE: RETURN,-1006 ; wrong type definition
END
    if (N_ELEMENTS(Mask) EQ 0) THEN     Mask=BYTE(0) ELSE    Mask=BYTE(Mask) 
    Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
    type_TMP=BYTARR(128)
    type_TMP(*)=0
    if (N_ELEMENTS(type) GT 0) THEN if (STRLEN(type) GT 0) THEN type_TMP(0:STRLEN(type)-1)=BYTE(type)
    command_TMP=BYTARR(128)
    command_TMP(*)=0
    if (N_ELEMENTS(command) GT 0) THEN if (STRLEN(command) GT 0) THEN command_TMP(0:STRLEN(command)-1)=BYTE(command)
IF (N_PARAMS() LT 6) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

MEM_TYP=[-1,1,2,4,4,8,8,-1,-1,16,-1,-1,-1,-1,8,-1]
 PIOWRITEVECTOBJECTMASK=call_external(PIOLibIDLSO,'piowritevectobjectmask_tempoidl', $
        data, $
        LONG64(N_ELEMENTS(data)*MEM_TYP(size(data,/type))), $
        Mask, $
        LONG64(N_ELEMENTS(Mask)*MEM_TYP(size(Mask,/type))), $
        Objectname_TMP, $
        type_TMP, $
        command_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOWRITEVECTOBJECTMASK
END
