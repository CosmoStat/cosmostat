
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetCRPixGrp

FUNCTION PIOGETCRPIXGRP,CRPix1,CRPix2,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(CRPix1) EQ 0) THEN     CRPix1=FLOAT(0) ELSE    CRPix1=FLOAT(CRPix1)
     if (N_ELEMENTS(CRPix2) EQ 0) THEN     CRPix2=FLOAT(0) ELSE    CRPix2=FLOAT(CRPix2)
 
 PIOGETCRPIXGRP=call_external(PIOLibIDLSO,'piogetcrpixgrp_tempoidl', $
        CRPix1, $
        CRPix2, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETCRPIXGRP
END
