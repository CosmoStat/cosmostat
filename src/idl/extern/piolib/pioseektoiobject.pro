
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSeekTOIObject

FUNCTION PIOSEEKTOIOBJECT,SampleNumber,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(SampleNumber) EQ 0) THEN     SampleNumber=LONG64(0) ELSE    SampleNumber=LONG64(SampleNumber)
 
 PIOSEEKTOIOBJECT=call_external(PIOLibIDLSO,'pioseektoiobject_tempoidl', $
        SampleNumber, $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOSEEKTOIOBJECT
END
