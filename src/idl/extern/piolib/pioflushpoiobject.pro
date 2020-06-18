
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFlushPOIObject

FUNCTION PIOFLUSHPOIOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIOFLUSHPOIOBJECT=call_external(PIOLibIDLSO,'pioflushpoiobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOFLUSHPOIOBJECT
END
