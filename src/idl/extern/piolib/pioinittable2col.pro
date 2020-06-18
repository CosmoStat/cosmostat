FUNCTION PIOINITTABLE2COL,NBPOINT,XUNIT,YUNIT,XNAME,YNAME,COMMENT, $
  XTYPE,YTYPE

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

TABLEUNIT=LONG64(0)

PIOERRMESS,call_external(PIOLibIDLSO,'pioinittable2col_tempoidl', $
              LONG64(NBPOINT),PIOSTRING2C(XUNIT),PIOSTRING2C(YUNIT), $
              PIOSTRING2C(XNAME),PIOSTRING2C(YNAME),PIOSTRING2C(COMMENT), $
              PIOSTRING2C(XTYPE),PIOSTRING2C(YTYPE),TABLEUNIT,/L64_VALUE)

CASE (STRING(XTYPE)) OF
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
        ELSE: return,-1
    END

CASE (STRING(YTYPE)) OF
        "PIOLONG": BEGIN
            Y=LONG64(0)
            DATAY = LON64ARR(NBPOINT)
            END
        "PIOBYTE": BEGIN
            Y=BYTE(0)
            DATAY = BYTARR(NBPOINT)
            END
        "PIOSHORT": BEGIN
            Y=0
            DATAY = INTARR(NBPOINT)
            END
        "PIOFLOAT": BEGIN
            Y=FLOAT(0)
            DATAY = FLTARR(NBPOINT)
            END
        "PIODOUBLE": BEGIN
            Y=DOUBLE(0)
            DATAY = DBLARR(NBPOINT)
            END
        "PIOCOMPLEX": BEGIN
            Y=COMPLEX(0,0)
            DATAY = COMPLEXARR(NBPOINT)
            END
        "PIODBLCOMPLEX": BEGIN
            Y=DCOMPLEX(0,0)
            DATAY = DCOMPLEXARR(NBPOINT)
            END
        ELSE: return,-1
    END


TABLE ={NBPOINT:NBPOINT, $
        XUNIT:XUNIT,YUNIT:YUNIT, $
        XNAME:XNAME,YNAME:YNAME, $
        COMMENT:COMMENT, XTYPE:XTYPE, $
        YTYPE:YTYPE,TABLEUNIT:TABLEUNIT, $
        X:DATAX,DX:DATAX,Y:DATAY,DY:DATAY}

; Free memory allocated for table

PIOERRMESS,call_external(PIOLibIDLSO,'piodeletetable2col_tempoidl', $
                        LONG64(TABLEUNIT))

RETURN,TABLE
END
