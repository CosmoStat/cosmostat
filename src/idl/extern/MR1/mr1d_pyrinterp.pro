;+
; NAME:
;        MR1D_PYRINTERP
;
; PURPOSE:
;	Interpolate each scale of a 1D multiresolution transform to the size
;	of the first scale
;
; CALLING:
;
;      MR1D_PYRINTERP, MR_Trans
;
; INPUTS-OUTPUTS:
;     MR_Trans -- 2D IDL array: Pyramidal multiresolution transform
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-

pro MR1D_PYRINTERP, MR_Trans

if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: MR1D_PYRINTERP, MR_Trans'
        goto, DONE
        end

Np = (size(MR_Trans))[1]
Nscale = (size(MR_Trans))[2]
TabNp = lonarr(Nscale)
TabNp[0] = Np
for j=1,Nscale-1 do TabNp[j] = (TabNp[j-1]+1)/2

for j=Nscale-1,0,-1 do BEGIN
  rec = MR_Trans[0:TabNp[j]-1, j]
;  for i=j-1,0,-1 do rec = congrid (rec, TabNp[i], /CUBIC)
  rec = congrid (rec, TabNp[0], /CUBIC)
  MR_Trans[*,j] = rec
  END

DONE:

END
