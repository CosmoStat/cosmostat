
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSetPOIObject

FUNCTION PIOSETPOIOBJECT,Tab,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
TYPE=BYTARR(128)
NBDATA= call_external(PIOLibIDLSO,'pioreadgetobjtype_idl', $
        MyObject, TYPE,/L64_VALUE) 
 CASE (STRING(TYPE)) OF
"PIOLONG": Tab = LONG64(Tab)
"PIOFLAG": Tab = BYTE(Tab)
"PIOBYTE": Tab = BYTE(Tab)
"PIOINT": Tab = LONG(Tab)
"PIOSTRING":BEGIN
    Tab=PIOSTRING2C(Tab)
     NBDATA=N_ELEMENTS(Tab)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(Tab(I))-1)=BYTE(Tab(I))
   Tab=RES
END 
"PIOSHORT": Tab = FIX(Tab)
 "PIOFLOAT": Tab = FLOAT(Tab)
 "PIODOUBLE": Tab = DOUBLE(Tab)
 "PIOCOMPLEX": Tab = COMPLEX(Tab)
 "PIODBLCOMPLEX": Tab = DBLCOMPLEX(Tab)
ELSE: RETURN,-1006 ; wrong type definition
END

MEM_TYP=[-1,1,2,4,4,8,8,-1,-1,16,-1,-1,-1,-1,8,-1]
 PIOSETPOIOBJECT=call_external(PIOLibIDLSO,'piosetpoiobject_tempoidl', $
        Tab, $
        LONG64(N_ELEMENTS(Tab)*MEM_TYP(size(Tab,/type))), $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOSETPOIOBJECT
END
