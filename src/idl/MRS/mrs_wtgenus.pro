;+
; NAME:
;        mrs_wtgenus
;
; PURPOSE:
;	 Computes the Genus curve for each scale of the Wavelet Transform.
;
; CALLING:
;
;     mrs_wtgenus, Image, Genus, NbrScale=NbrScale, Pyr=Pyr, Trans=Trans, MFSTEP=MFSTEP        
;    
; INPUT:
;     Imag -- IDL array of healpix map: Input image be filtered 
;
; OUTPUT:
;     Genus -- 3D IDL array [*,2,NbrScale-1]: Genus curve at each scale
;                            Genus[*,0,j] = Nu value of the Genus curve at scale j
;                            Genus[*,1,j] = Genus values at scale j   
;
; INPUT KEYWORDS:
;      NbrScale -- int: Number of scales (default is 4)
;      Pyr      -- int: if set, use of the pyramidal WT instead of the undecimated one.
;      MFSTEP   -- float: step parameter for the evaluation of the Genus. Default is 0.005
;
; OUTPUT KEYWORDS:
;      Trans -- IDL structure: Wavelet transform of the data
;
; EXTERNAL CALLS:
;       mrs_wttrans
;       mrs_pwttrans
;
; EXAMPLE:
;     npix = 196608l
;     n = randomn(seed, npix)
;     mrs_wtgenus, n, Genus, NbrScale=5, /pyr 
;     ; Plot the genus curve at each scale j
;     for j=0, 3 do plot, Genus[*,0,j], Genus[*,1,j]  
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2004
;	February, 2005 File creation
;-
;=================================================================================================================

function get_genus_val, Ima, Threshold, NuXAbs, NbrIsol=NbrIsol, NbrHole=NbrHole, GS=GS

ThresIma = Ima
ThresIma[*]= 0
ind = where( Ima GE Threshold, c)
if c GT 0 then begin
   ThresIma[ind] = 1
   l = long(label_region(ThresIma, /ALL_NEIGHBORS, /ULONG))
   NbrIsol = long(max(l))
end else NbrIsol = 0

if c GT 0 then  begin
   ThresIma[*]= 1
   ThresIma[ind] = 0
   l =  long(label_region(ThresIma, /ALL_NEIGHBORS, /ULONG))
   NbrHole =  long(max(l))
end else NbrHole = 0

;Genus = long(NbrIsol) - long(NbrHole) + 1l

Genus = 0l
if NuXAbs GT 0 then  begin
   Genus = NbrHole - NbrIsol + 1l
   GS =  NbrIsol 
   end else  begin
   Genus  = NbrIsol - NbrHole + 1l
   GS  =  NbrHole
end 

; print, Threshold, NbrIsol, NbrHole, Genus

return, Genus
end

;=====================================================================================================================

function xerf, X
    Z = X;
    if X GT 0.5 then Z = 1 - X
    if Z LT 1e-10 then  Val_Return  = 0. $
    else begin
       Z = sqrt (-2d * alog (Z))
       A0 = 2.30753d
       A1 = 0.27061d
       B1 = 0.99229d
       B2 = 0.04481d
       Val_Return = Z - (A0 + A1 * Z) / (1 + ((B1 + B2 * Z) * Z))
       if X GT 0.5 then Val_Return = - Val_Return
    end
    return, -Val_Return
end

;======================================================================================================================

pro tg, GenFun

np=256
n = randomn(seed, np,np)
g = mygauss(np,np,7)
gn = float( dfti(dft(n)*dft(g)))
gn = gn / sigma(gn)
info, gn
genus2d, gn, GenFun
plot, GenFun(*,0), GenFun(*,1)

end

;=======================================================================================================================

pro genus5sig_2d, Ima, Ima1, GenFun, GenFun1, SigmaStep=SigmaStep, NSigma=NSigma, nstruc=nstruc

if not keyword_set(SigmaStep) then SigmaStep = 0.1
if not keyword_set(NSigma) then NSigma = 5.

 
Np = long(2.*NSigma/SigmaStep+1)
if keyword_set(nstruc) then GenFun = fltarr(Np, 3) $
else GenFun = fltarr(Np, 2)
Npix = long(N_ELEMENTS(Ima))
if keyword_set(nstruc) then GenFun1 = fltarr(Np, 3) $
else GenFun1 = fltarr(Np, 2)
sdens = Ima
inds = sort(Ima)
sdens[*] = Ima[inds] 
p=0
sig=sigma(ima)
for vf=-NSigma,NSigma,SigmaStep  do begin
    iro=long(vf*Npix)
    G=get_genus_val(Ima, vf*sig, vf, GS=GS)
    G1=get_genus_val(Ima1, vf*sig, vf, GS=GS1)
    GenFun(p, 0) = vf;
    GenFun(p, 1) = G
    GenFun1(p, 0) = vf;
    GenFun1(p, 1) = G1
    if keyword_set(nstruc) then GenFun(p, 2) = GS
    if keyword_set(nstruc) then GenFun1(p, 2) = GS1
    p = p + 1
end
GenFun=GenFun(0:p-1,*)
GenFun1=GenFun1(0:p-1,*)

GenFun[*,1] = 2.*(1.- GenFun[*,1] )
GenFun1[*,1] = 2.*(1.- GenFun1[*,1] )

end

;G = 1-M/2
;G-1=-M/2
;M= 2*(1-G)

;===========================================================================================================

pro genus2d, Ima, GenFun, MFSTEP=MFSTEP, nstruc=nstruc

if not keyword_set(MFSTEP) then MFSTEP = 0.005
VFSTEP=MFSTEP
Np = long(1./VFSTEP)
if keyword_set(nstruc) then GenFun = fltarr(Np, 3) $
else GenFun = fltarr(Np, 2)
Npix = long(N_ELEMENTS(Ima))

sdens = Ima
inds = sort(Ima)
sdens[*] = Ima[inds] 
p=0
for vf=VFSTEP,1.-VFSTEP,VFSTEP  do begin
    iro=long(vf*Npix)
    ro=sdens(iro)
    NuXAbs = xerf(vf)
    G=get_genus_val(Ima, ro,NuXAbs, GS=GS)
    GenFun(p, 0) = NuXAbs;
    GenFun(p, 1) = G
    if keyword_set(nstruc) then GenFun(p, 2) = GS
    p = p + 1
end
GenFun = GenFun(0:p-1,*)
end

;==============================================================================================================

pro healpixgenus2d, Ima, GenFun, MFSTEP=MFSTEP

if not keyword_set(MFSTEP) then MFSTEP = 0.005
VFSTEP=MFSTEP
Np = long(1./VFSTEP)
GenFun = fltarr(Np, 2)
TabIsol = lonarr(12)
TabHole = lonarr(12)
Npix = long(N_ELEMENTS(Ima))

sdens = Ima
inds = sort(Ima)
sdens[*] = Ima[inds] 
p = 0l
get_all_faces, Ima, CubeFace

for vf=VFSTEP,1.-VFSTEP,VFSTEP  do begin
    iro=long(vf*Npix)
    ro=sdens(iro)
    NuXAbs = xerf(vf) 
    GenFun[p,0] = NuXAbs 
    
    for c=0, 11 do begin
       Face = CubeFace[*,*,c]
       G=get_genus_val(Face, ro, NuXAbs, NbrIsol=NbrIsol, NbrHole=NbrHole)
       TabIsol[c] = NbrIsol
       TabHole[c]  = NbrHole
    end
    GenFun[p,1] = total(TabIsol)  - total(TabHole)  + 1l
    p = p + 1
end
GenFun = GenFun[0:p-1,*]
; plot, GenFun(*,0), GenFun(*,1)
end

;========================================================================================================================


pro mrs_wtgenus, Imag, TabGen, NbrScale=NbrScale, Pyr=Pyr, MFSTEP=MFSTEP, Trans=Trans

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_wtgenus, Imag, TabGen, NbrScale=NbrScale, Pyr=Pyr, MFSTEP=MFSTEP, Trans=Trans'
        goto, DONE
        end
	
;if keyword_set(Mask) then begin
;   vs = size(Mask)
;   Nm = vs[1]
;   if Nm NE npix then begin
;       print, 'Error: the mask must have the same dimension as the input map'
;       goto , DONE
;   end
;end


if not keyword_set(NbrScale) then NbrScale=4
if not keyword_set(MFSTEP) then MFSTEP = 0.005

if not keyword_set(pyr) then mrs_wttrans, Imag, Trans, NbrScale=NbrScale $
else mrs_pwttrans, Imag, Trans, NbrScale=NbrScale

for j =0,NbrScale-2 do begin
    Scale = mrs_wtget(Trans,j,NormVal=NormVal)
    ; help, scale
    Scale = Scale / sigma(Scale)
    healpixgenus2d, Scale, GenFun, MFSTEP=MFSTEP   
    ; help, GenFun
    ; plot, GenFun(*,0), GenFun(*,1)
    if j eq 0 then begin
      vs = size(GenFun)
      Np = vs[1]
      TabGen = fltarr(Np,2,NbrScale-1)
      TabGen[*,*,j] = GenFun
    end  
    TabGen[*,*,j] = GenFun 
endfor

DONE:

END

