;+
; NAME:
;        MR1D_PYRMED
;
; PURPOSE:
;	Computes a multiresolution transform of a 1D Signal by the
;	min-max algorithm.
;
; CALLING:
;
;      MR1D_MINMAX, Signal, MinMaxTrans, Nscale
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want transform
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     MinMaxTrans -- 2D: multiresolution transform
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR1D_MINMAX, Signal, MinMaxTrans, Nscale

if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: MR1D_PYRMED, Signal,  MinMaxTrans, Nscale '
        goto, DONE
        end

Np = (size(Signal))[1]
MinMaxTrans = fltarr(Np, Nscale)

MinMaxTrans[*, 0] = Signal
MM = fltarr(Np)

Win=2
for j=0,Nscale-2 do $
BEGIN
    Scale = MinMaxTrans[*, j]
    i = long(0)
    while i LT Np-1 do $
    BEGIN
       voisin=Scale( max([i-Win,0]) : min([i+Win, Np-1]))
       MM[i] = min(voisin)
       i = i +1
    END
    MinMaxTrans[*, j] = MinMaxTrans[*, j]  - MM
    MinMaxTrans[*, j+1] = MM
    Win = 2*Win
END

DONE:

END
