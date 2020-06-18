;+
; NAME:
;        MR1D_PYRREC
;
; PURPOSE:
;	Reconstruct a signal from its pyramidal multiresolution transform
;
; CALLING:
;
;      MR1D_PYRREC, MR_Trans, Rec_Signal
;       
;
; INPUTS:
;     MR_Trans -- 2D IDL array: Pyramidal multiresolution transform
;
; OUTPUTS:
;     Rec_Signal -- 1D array: reconstructed signal
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-

pro MR1D_PYRREC, MR_Trans, Rec_Signal

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: MR1D_PYRREC, MR_Trans, Rec_Signal'
        goto, DONE
        end

Np = (size(MR_Trans))[1]
Nscale = (size(MR_Trans))[2]
TabNp = lonarr(Nscale)
TabNp(0) = Np

for j=1,Nscale-1 do TabNp[j] = (TabNp[j-1]+1)/2

j = Nscale-2
Rec_Signal = MR_Trans[0:TabNp[j+1]-1, j+1]

for j=Nscale-2,0,-1 do BEGIN
   Rec_Signal = congrid (Rec_Signal, TabNp[j], /CUBIC)
   Rec_Signal[*] = Rec_Signal[*] + MR_Trans[0:TabNp[j]-1, j]
   END

DONE:

END
