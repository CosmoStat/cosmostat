
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOKeywordListGrp

FUNCTION PIOKEYWORDLISTGRP,Keyword,Type,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Keyword=long64(0)
    Type=long64(0)

 PIOKEYWORDLISTGRP=call_external(PIOLibIDLSO,'piokeywordlistgrp_tempoidl', $
        Keyword, $
        Type, $
        MyGroup, $
               /L64_VALUE) 
  IF (PIOKEYWORDLISTGRP GT 0) THEN BEGIN
     TP=Keyword
     TMPNAME=bytarr(128)
     Keyword=STRARR(PIOKEYWORDLISTGRP)
        FOR I=0L,PIOKEYWORDLISTGRP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            Keyword(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (PIOKEYWORDLISTGRP GT 0) THEN BEGIN
     TP=Type
     TMPNAME=bytarr(128)
     Type=STRARR(PIOKEYWORDLISTGRP)
        FOR I=0L,PIOKEYWORDLISTGRP-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            Type(I)=strmid(string(TMPNAME),0,A)
            END
    END

 RETURN,PIOKEYWORDLISTGRP
END
