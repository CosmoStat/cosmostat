
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFlushTOIObject

FUNCTION PIOFLUSHTOIOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIOFLUSHTOIOBJECT=call_external(PIOLibIDLSO,'pioflushtoiobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOFLUSHTOIOBJECT
END
