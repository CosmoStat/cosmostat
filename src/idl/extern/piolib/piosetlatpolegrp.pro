
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSetLatpoleGrp

FUNCTION PIOSETLATPOLEGRP,LatPole,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(LatPole) EQ 0) THEN     LatPole=DOUBLE(0) ELSE    LatPole=DOUBLE(LatPole)
 
 PIOSETLATPOLEGRP=call_external(PIOLibIDLSO,'piosetlatpolegrp_tempoidl', $
        LatPole, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOSETLATPOLEGRP
END
