
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOFlushTIMESTAMPObject

FUNCTION PIOFLUSHTIMESTAMPOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIOFLUSHTIMESTAMPOBJECT=call_external(PIOLibIDLSO,'pioflushtimestampobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOFLUSHTIMESTAMPOBJECT
END
