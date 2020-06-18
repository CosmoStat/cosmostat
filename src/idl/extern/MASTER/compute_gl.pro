

PRO compute_gl, mll, mbb, matp, matq, lmax, nbins, mask=mask, filemask=filemask, filebinmax=filebinmax, filebinmin=filebinmin, deltl=deltl, saveout=saveout


;#########################################################################
;#  input:                                                              ;#
;#    DELTL is the number of ells per bin.                              ;#
;#        used only if FILEBINMAX and FILEBINMIN are not set            ;#
;#    MASK is the healpix map of the mask (or the pixel weighting)      ;#
;#    lmax: maximum value of ell for the coupling matrix computation    ;# 
;#  output:                                                             ;#
;#    mll is the coupling matrix                                        ;#
;#    mbb is the inverse of the binned coupling matrix                  ;#
;#    matp and matq are the "binning" matrices such that                ;#
;#                     cb = matp#cl                                     ;#
;#                     cl = matq#cb (assuming flat-band power spectrum) ;#
;#  G.P. 2004                                                           ;#
;#########################################################################


IF (KEYWORD_SET(filemask) EQ 1) THEN BEGIN
    read_fits_map, filemask, mapmask
    PRINT, 'file ', filemask, ' is read'
ENDIF
IF (KEYWORD_SET(mask) EQ 1) THEN mapmask=mask

nside = LONG(sqrt(N_ELEMENTS(mapmask)/12))
print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print, 'nside=',nside
print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

;lecture des fichiers concernant le binning
IF KEYWORD_SET(filebinmax) EQ 1 THEN BEGIN
    close, 1
    nlmax = dblarr(nbins)
    openr, 1, filebinmax
    readf, 1, nlmax
    close, 1
    IF KEYWORD_SET(filebinmin) EQ 1 THEN BEGIN
        close, 1
        nlmax = dblarr(nbins)
        openr, 1, filebinmin
        readf, 1, nlmin
        close, 1
    ENDIF ELSE BEGIN
        nlmin = nlmax+1
        nlmin = [0, nlmin(0:nbins-2)]
    ENDELSE
ENDIF ELSE BEGIN
    IF KEYWORD_SET(filebinmin) EQ 1 THEN BEGIN
        close, 1
        nlmax = dblarr(nbins)
        openr, 1, filebinmin
        readf, 1, nlmin
        close, 1
        nlmax = nlmin-1
        nlmax = [nlmax(1:nbins-1), lmax]
    ENDIF ELSE BEGIN
        nlmax = lindgen(nbins)*deltl+deltl-1
        nlmin = lindgen(nbins)*deltl
    ENDELSE
ENDELSE

IF (nlmin(nbins-1) LE lmax) THEN BEGIN
    nlmax(nbins-1) = lmax
ENDIF ELSE BEGIN
    print, 'ERROR: check lmax'
    nlmin='NAN'
ENDELSE

;compute the power spectrum of the mask map
;master_anafast_h, mapmask, ellw, well, 3*nside-1, tmpdir='/tmp', diranafast='';/home/patanch/Projets/Archeops/Softs/'

well = mrs_powspec(mapmask,/nonorm)
ellw = findgen(n_elements(well))


;put 0 for l>(3*nside+1) 
;well2 = DBLARR(2l*(3*nside-1))
;ellw = findgen(2l*lmax+1)
;well2(0:(3*nside-1)) = well(0:(3*nside-1))
;well = well2(0:2*lmax)

;compute the coupling matrix
master_coupling_matrix,ellw,well,lmax,deltal,mll,ell,ellbins,matp,matq,bintabmax=nlmax,bintabmin=nlmin, /csbin


kmat=master_compute_cormat(mll,matp,matq)
kmatinv=invert(kmat)

Mbb = kmatinv#matp

IF KEYWORD_SET(saveout) EQ 1 THEN save, file=saveout, mll, Mbb, matp, matq


RETURN

END
