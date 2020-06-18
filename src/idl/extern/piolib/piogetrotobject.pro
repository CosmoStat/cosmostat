
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetRotObject

FUNCTION PIOGETROTOBJECT,Rot[9],Objectname,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 CASE (STRING(TYPE)) OF
"PIOLONG": Rot[9] = LONG64(Rot[9])
"PIOFLAG": Rot[9] = BYTE(Rot[9])
"PIOBYTE": Rot[9] = BYTE(Rot[9])
"PIOINT": Rot[9] = LONG(Rot[9])
"PIOSTRING":BEGIN
    Rot[9]=PIOSTRING2C(Rot[9])
     NBDATA=N_ELEMENTS(Rot[9])
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(Rot[9](I))-1)=BYTE(Rot[9](I))
   Rot[9]=RES
END 
"PIOSHORT": Rot[9] = FIX(Rot[9])
 "PIOFLOAT": Rot[9] = FLOAT(Rot[9])
 "PIODOUBLE": Rot[9] = DOUBLE(Rot[9])
 "PIOCOMPLEX": Rot[9] = COMPLEX(Rot[9])
 "PIODBLCOMPLEX": Rot[9] = DBLCOMPLEX(Rot[9])
ELSE: RETURN,-1006 ; wrong type definition
END
    Objectname_TMP=BYTARR(128)
    Objectname_TMP(*)=0
    if (N_ELEMENTS(Objectname) GT 0) THEN if (STRLEN(Objectname) GT 0) THEN Objectname_TMP(0:STRLEN(Objectname)-1)=BYTE(Objectname)
IF (N_PARAMS() LT 3) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOGETROTOBJECT=call_external(PIOLibIDLSO,'piogetrotobject_tempoidl', $
        Rot[9], $
        Objectname_TMP, $
        MyGroup, $
               /L64_VALUE) 
  Objectname=STRING(Objectname_TMP)

 RETURN,PIOGETROTOBJECT
END
