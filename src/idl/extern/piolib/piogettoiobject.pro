
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTOIObject

FUNCTION PIOGETTOIOBJECT,data,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
TYPE=BYTARR(128)
NBDATA= call_external(PIOLibIDLSO,'pioreadgetobjtype_idl', $
        MyObject, TYPE,/L64_VALUE) 
CASE (STRING(TYPE)) OF
  "PIOFLAG": DATA = BYTARR(NBDATA)
  "PIOLONG": DATA = LON64ARR(NBDATA)
  "PIOBYTE": DATA = BYTARR(NBDATA)
  "PIOINT": DATA = LONARR(NBDATA)
  "PIOSTRING": DATA = BYTARR(128L*NBDATA)
  "PIOSHORT": DATA = INTARR(NBDATA)
  "PIOFLOAT": DATA = FLTARR(NBDATA)
  "PIODOUBLE": DATA = DBLARR(NBDATA)
  "PIOCOMPLEX": DATA = COMPLEXARR(NBDATA)
  "PIODBLCOMPLEX": DATA = DCOMPLEXARR(NBDATA)
END
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

 PIOGETTOIOBJECT=call_external(PIOLibIDLSO,'piogettoiobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 
IF (PIOGETTOIOBJECT NE 0) THEN $
A=call_external(PIOLibIDLSO,'piomemcpy_nofreeidl',PIOGETTOIOBJECT,DATA, LONG64(NBDATA*PIOTYPESIZE(STRING(TYPE))))

 RETURN,PIOGETTOIOBJECT
END
