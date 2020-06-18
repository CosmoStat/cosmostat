;+
; NAME: 
;       read_cambcl
;
; PURPOSE: 
;      read the Cl from a file in Camb format. 
;
; CALLING:
;      Cl = read_cambcl(FileName, ell=ell, LL1=LL1, nomonodip=nomonodip)
;
; INPUT:
;       FileName: File Name 
;
; INPUT KEYWORDS:
;        LL1: scalar --  the returned Cl are multiplied by l(l+1)/2PI, i.e. as in standard CAMB output
;        nomonodip: scalar -- remove the monopole and the dipole, i.e. as in standard CAMB output
;                             otherwise add the C_0 and C_1 values.
;
; OUTPUT KEYWORDS:
;      Ell: IDL 1D array:  l values related to the CL
;
; OUTPUTS:
;           Cl: IDL 1D array -- Cl values
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;           Cl = read_cambcl(FileName)
;            
;           Cl = read_cambcl(FileName, /LL1, /nomonodip)
;           ; Read the ouput CAMB exactly as it is (without adding a monopole and a dipole and without renormalization with l(l+1)/2PI).
;
; HISTORY:
;	Written: Jean-Luc Starck 2013
;    January 2014 - Update the default Cl to Planck-PR1 official Cl best fit (/planck_lowl/base_planck_lowl.bestfit_cl  in ESA Planck PAIO)
;-
;==================================================
;   FN = '/Users/starck/Main/ISAP/param/CLASScmb/PlanckParams/Planck/pr1classcl_lensed.dat'
;    FN = '/Users/starck/Main/ISAP/param/cmb/PlanckParams/Planck/base_planck_lowl/base/planck_lowl/base_planck_lowl.bestfit_cl'
;    z = read_cambcl(FN)
;    z = z [0:3000]

;    FN1 = '/Users/starck/Main/ISAP/param/cmb/PlanckParams/PlanckWP/base_planck_lowl_lowLike/base/planck_lowl_lowLike/base_planck_lowl_lowLike.bestfit_cl'
;    z1 = read_cambcl(FN1)
;    z1 = z1 [0:3000]
;    plot, (z-z1)/z * 100.

;  j  = rim('~/CambCl/best_fits_pr1.fits')
;  x = getcmb(cl=T0, /no, /OLD)
;  x = getcmb(cl=T1, /no)
;  plot, (T0-T1)/T0 * 100.

;  plot, (T0-T1)/T0 * 100.
;  plot, (T1-z)/T1 * 100.
;  plot, (T1-j)/T1 * 100.
;  plot, (z-j)/z * 100.

; p = read_cambcl('/Users/starck/CambCl/P_BF/_lensedCls.dat')
;  plot, (z-p)/z * 100.
;  plot, (T1-p)/T1 * 100.

function read_cambcl, FileName, LL1=LL1, ell=ell, nomonodip=nomonodip, DATA_START=DATA_START
if not keyword_set(DATA_START) then DATA_START = 0
; print, FileName
a = read_ascii(FileName, DATA_START=DATA_START)
if a.field1[0,0] NE 2 then print, 'WARNING: the Cl vector does not start at l=2, something wrong happened … '
a=a.field1
ell = long(reform(a[0,*]))
Cl = double(reform(a[1, *]))
L1 = ell*(ell+1) / (2. * !DPI)
if not keyword_set(ll1) then cl = cl / L1
if not keyword_set(nomonodip) then begin
Cl = [0,0, Cl]
ell = [0,1, ell]
end
return, Cl
end

;===================================================
; #   Expected input format    L    TT             TE             EE             BB             PP
function read_cambclp, FileName, LL1=LL1, ell=ell, nomonodip=nomonodip, DATA_START=DATA_START
if not keyword_set(DATA_START) then DATA_START = 0
; print, FileName
a = read_ascii(FileName, DATA_START=DATA_START)

if a.field1[0,0] NE 2 then print, 'WARNING: the Cl vector does not start at l=2, something wrong happened … '
a=a.field1
ell = long(reform(a[0,*]))
Cl_TT = double(reform(a[1, *]))
Cl_TE = double(reform(a[2, *]))
Cl_EE = double(reform(a[3, *]))
Cl_BB = double(reform(a[4, *]))
Cl_PP = double(reform(a[5, *]))
 
L1 = ell*(ell+1) / (2. * !DPI)
if not keyword_set(ll1) then begin
   Cl_TT = Cl_TT / L1
   Cl_TE = Cl_TE / L1
   Cl_EE = Cl_EE / L1
   Cl_BB = Cl_BB / L1
   Cl_PP = Cl_PP / L1
end

if not keyword_set(nomonodip) then begin
Cl_TT = [0,0, Cl_TT]
Cl_TE = [0,0, Cl_TE]
Cl_EE = [0,0, Cl_EE]
Cl_BB = [0,0, Cl_BB]
Cl_PP = [0,0, Cl_PP]
ell = [0,1, ell]
; Healpix format ;  (lmax+1,6) given in the sequence T E B TxE TxB ExB
;                    (lmax+1,4) given in the sequence T E B TxE,
;                     for temperature and polarisation
Cl = [[Cl_TT], [Cl_EE], [Cl_BB], [Cl_TE]]
end
return, Cl
end

