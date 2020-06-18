;+
; NAME:
;        MR1D_PAVEMED
;
; PURPOSE:
;	Computes a multiresolution transform of a 1D Signal by the
;	multiresolution median algorithm.
;
; CALLING:
;
;      MR1D_PAVEMED, Signal, MedTrans, Nscale
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want transform
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     MedTrans -- 2D IDL array: multiresolution transform
;
; INPUT KEYWORDS:
;     cont -- scalar: if set, consider the data is constant after
;                     the borders
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR1D_PAVEMED, Signal, MedTrans, Nscale, cont=cont

if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: MR1D_PAVEMED, Signal,  MedTrans, Nscale, cont=cont'
        goto, DONE
        end

Np = (size(Signal))[1]
MedTrans = fltarr(Np, Nscale)

MedTrans(*, 0) = Signal
Win = fltarr(5)
med =  fltarr(Np)

cont=1
Win=2
for j=0,Nscale-2 do $
BEGIN
    Scale = MedTrans[* , j]
    i = long(0)
    while i LT Np do $
    BEGIN
       if keyword_set(cont) and ((i-Win LT 0) or (i+ 2 GT Np-1)) then $
       BEGIN
          if i - Win LT 0 then voisin = Scale(0:2*Win+1) $
          else if i + Win GT Np-1 then voisin = Scale(Np-1 - 2*Win+1: Np-1) 
       END ELSE voisin=Scale( max([i-Win,0]) : min([i+Win, Np-1]))
       med[i] = voisin( (sort(voisin))[ n_elements(voisin)/2 ])
       i = i +1
    END
    MedTrans[*, j+1] = med
    MedTrans[*, j] = Scale[*] - med[*]
    Win = 2*Win
END
DONE:

END
