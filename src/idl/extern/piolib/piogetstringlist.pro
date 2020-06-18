FUNCTION PIOGETSTRINGLIST,VAL,NAME,ATTRIB,IMOUNIT

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

TEMP=bytarr(128)

TPLIST = LONG64(0)
SZ=LONG64(0)
SZ= call_external(PIOLibIDLSO,'piogetstringlist_tempoidl', $
                      TPLIST,PIOSTRING2C(NAME),PIOSTRING2C(ATTRIB), $
                    LONG64(IMOUNIT),/L64_VALUE)

IF (sz gt 0) THEN BEGIN
    VAL = STRARR(SZ)
    FOR I=0,SZ-1 DO BEGIN
        A=call_external(PIOLibIDLSO,'piogetstringidl',TPLIST+I*128,TEMP)
        VAL(I)=strmid(STRING(TEMP),0,A)
    END
END
RETURN,sz
END
