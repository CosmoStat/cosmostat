FUNCTION PIOGETTABLE,NAME,IMOUNIT

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

TEMPNAME=BYTARR(128)
TEMPNAME1=TEMPNAME
TEMPNAME2=TEMPNAME
TEMPNAME3=TEMPNAME
SZ0=LONG64(0)
SZ1=LONG64(0)
SZ2=LONG64(0)
SZ3=LONG64(0)

ERR = LONG64(0)
NBPOINT= LONG64(0)
TABLEUNIT= LONG64(0)

ERR = call_external(PIOLibIDLSO,'piogettable_tempoidl', $
              LONG64(NBPOINT),TEMPNAME,TEMPNAME1,TEMPNAME2,TEMPNAME3, $
              TABLEUNIT, $
              PIOSTRING2C(NAME),LONG64(IMOUNIT), $
              SZ0,SZ1,SZ2,SZ3,/L64_VALUE)
PIOERRMESS,ERR

TYPE =STRMID(STRING(TEMPNAME3),0,SZ2)
CASE (STRING(TYPE)) OF
        "PIOLONG": BEGIN
            X=LONG64(0)
            DATAX = LON64ARR(NBPOINT)
            END
        "PIOBYTE": BEGIN
            X=BYTE(0)
            DATAX = BYTARR(NBPOINT)
            END
        "PIOSHORT": BEGIN
            X=0
            DATAX = INTARR(NBPOINT)
            END
        "PIOFLOAT": BEGIN
            X=FLOAT(0)
            DATAX = FLTARR(NBPOINT)
            END
        "PIODOUBLE": BEGIN
            X=DOUBLE(0)
            DATAX = DBLARR(NBPOINT)
            END
        "PIOCOMPLEX": BEGIN
            X=COMPLEX(0,0)
            DATAX = COMPLEXARR(NBPOINT)
            END
        "PIODBLCOMPLEX": BEGIN
            X=DCOMPLEX(0,0)
            DATAX = DCOMPLEXARR(NBPOINT)
            END
        ELSE: return,0
    END


TABLE ={NBPOINT:NBPOINT,UNIT:'',NAME:'', $
        COMMENT:'', TYPE:'',TABLEUNIT:LONG64(0), $
        X:DATAX,DX:DATAX}

TABLE.UNIT = STRMID(STRING(TEMPNAME),0,SZ0)
TABLE.COMMENT =STRMID(STRING(TEMPNAME1),0,SZ1)
TABLE.NAME =STRMID(STRING(TEMPNAME2),0,SZ3)
TABLE.TYPE =TYPE

for i=0,NBPOINT-1 do begin
    DX=X
    err = call_external(PIOLibIDLSO,'pioextracttable_tempoidl', $
                        X,DX,PIOSTRING2C(TABLE.TYPE), $
                        LONG64(i), LONG64(TABLEUNIT))
    PIOERRMESS,ERR
    TABLE.X(I)=X
    TABLE.DX(I)=DX
end

; Free memory allocated for table

PIOERRMESS,call_external(PIOLibIDLSO,'piodeletetable_tempoidl', $
                        LONG64(TABLEUNIT))
RETURN,TABLE
END
