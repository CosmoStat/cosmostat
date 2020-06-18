;+
; NAME:
;        MR1D_PYRMED
;
; PURPOSE:
;	Computes a multiresolution transform of a 1D Signal by the
;	pyramidal median algorithm.
;
; CALLING:
;
;      MR1D_PYRMED, Signal, MedTrans, Nscale, TabNp=TabNp, interp=interp,  
;                   cont=cont
;       
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want transform
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     MedTrans -- 2D: multiresolution transform
;
; INPUT KEYWORDS:
;     interp -- scalar: if set, each scale is interpolated to the size of the 
;                       input signal.
;     cont -- scalar: if set, consider the data is constant after
;                     the borders
;
; OUTPUT KEYWORDS:
;    TabNp -- 1D array: size of each scale
;
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR1D_PYRMED, Signal, MedTrans, Nscale, TabNp=TabNp, interp=interp,cont=cont

if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: MR1D_PYRMED, Signal,  MedTrans, Nscale '
        print, '                 TabNp=TabNp, interp=interp, cont=cont'
        goto, DONE
        end

Np = (size(Signal))[1]
MedTrans = fltarr(Np, Nscale)
TabNp = lonarr(Nscale)
TabNp[0] = Np

MedTrans[*, 0] = Signal
Win = fltarr(5)

for j=0,Nscale-2 do $
BEGIN
    Scale = MedTrans[0:TabNp[j]-1 , j]
    med =  fltarr(TabNp[j])
    i = long(0)
    while i LT TabNp(j)-1 do $
    BEGIN
       voisin=Scale[ max([i-2,0]) : min([i+2, TabNp(j)-1])]
       if keyword_set(cont) and ((i-2 LT 0) or (i+ 2 GT TabNp(j)-1)) then $
       BEGIN
          if i - 2 LT 0 then Valcont = median(Signal[0:15]) $
          else if i + 2 GT TabNp[j]-1 then Valcont = median(Scale[TabNp[j]-16:TabNp[j]-1])
          Win[*] = Valcont
          Win[0: (size(voisin))[1]-1] = voisin
       END ELSE med[i] = voisin[ (sort(voisin))[ n_elements(voisin)/2 ]]
       i = i +1
    END
    TabNp[j+1] = (TabNp[j]+1)/2
    reb = congrid(med, TabNp[j+1] )
    MedTrans[0:TabNp[j+1]-1, j+1] = reb[*]

    med = congrid (reb, TabNp[j], /CUBIC)
    MedTrans[0:TabNp[j]-1, j] = Scale[*] - med[*]
END

  if keyword_set(interp) then mr1d_pyrinterp, MedTrans

DONE:

END
