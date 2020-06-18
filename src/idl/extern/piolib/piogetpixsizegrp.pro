
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetPixSizeGrp

FUNCTION PIOGETPIXSIZEGRP,PixSize1,PixSize2,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(PixSize1) EQ 0) THEN     PixSize1=FLOAT(0) ELSE    PixSize1=FLOAT(PixSize1)
     if (N_ELEMENTS(PixSize2) EQ 0) THEN     PixSize2=FLOAT(0) ELSE    PixSize2=FLOAT(PixSize2)
 
 PIOGETPIXSIZEGRP=call_external(PIOLibIDLSO,'piogetpixsizegrp_tempoidl', $
        PixSize1, $
        PixSize2, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETPIXSIZEGRP
END
