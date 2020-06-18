
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSetTIMESTAMPObject

FUNCTION PIOSETTIMESTAMPOBJECT,TabIn,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
TYPE=BYTARR(128)
NBDATA= call_external(PIOLibIDLSO,'pioreadgetobjtype_idl', $
        MyObject, TYPE,/L64_VALUE) 
 CASE (STRING(TYPE)) OF
"PIOLONG": TabIn = LONG64(TabIn)
"PIOFLAG": TabIn = BYTE(TabIn)
"PIOBYTE": TabIn = BYTE(TabIn)
"PIOINT": TabIn = LONG(TabIn)
"PIOSTRING":BEGIN
    TabIn=PIOSTRING2C(TabIn)
     NBDATA=N_ELEMENTS(TabIn)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(TabIn(I))-1)=BYTE(TabIn(I))
   TabIn=RES
END 
"PIOSHORT": TabIn = FIX(TabIn)
 "PIOFLOAT": TabIn = FLOAT(TabIn)
 "PIODOUBLE": TabIn = DOUBLE(TabIn)
 "PIOCOMPLEX": TabIn = COMPLEX(TabIn)
 "PIODBLCOMPLEX": TabIn = DBLCOMPLEX(TabIn)
ELSE: RETURN,-1006 ; wrong type definition
END

MEM_TYP=[-1,1,2,4,4,8,8,-1,-1,16,-1,-1,-1,-1,8,-1]
 PIOSETTIMESTAMPOBJECT=call_external(PIOLibIDLSO,'piosettimestampobject_tempoidl', $
        TabIn, $
        LONG64(N_ELEMENTS(TabIn)*MEM_TYP(size(TabIn,/type))), $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOSETTIMESTAMPOBJECT
END
