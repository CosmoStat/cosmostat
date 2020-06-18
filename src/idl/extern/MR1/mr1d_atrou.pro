;+
; NAME:
;        MR1D_ATROU
;
; PURPOSE:
;	Computes a multiresolution transform of a 1D Signal by the
;	a-trous algorithm.
;
; CALLING:
;
;      MR1D_ATROU, Signal, WaveTrans, Nscale
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want transform
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     WaveTrans -- 2D: multiresolution transform
;
; KEYWORDS:
;     mirror -- scalar: if set, then a mirror is used at the border
;
;     linear -- scalar: if set, then a linear wavelet function is used.
;                        if not set, then B-spline wavelet function is used.
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR1D_ATROU, Signal, WaveTrans, Nscale, mirror=mirror, linear=linear

if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: MR1D_ATROU, Signal,  WaveTrans, Nscale, mirror=mirror, linear=linear'
        goto, DONE
        end

Np = (size(Signal))[1]
WaveTrans = fltarr(Np, Nscale)
Scale = fltarr(Np)

WaveTrans[*, 0] = Signal

for j=0,Nscale-2 do $
BEGIN
   step = 2^j
   i = long(0)
   while i LT Np do $
      BEGIN
         if keyword_set(mirror) then BEGIN
            im1 = i-step
            if im1 lt 0 then im1 = - im1
            im2 = i-2*step
            if im2 lt 0 then im2 = - im2
            ip1 = i+step
            if ip1 ge Np then ip1 = 2*Np-2-ip1
            ip2 = i+2*step
            if ip2 ge Np then ip2 = 2*Np-2-ip2
          END ELSE BEGIN
            im1 = max([i-step,0])
            im2 = max([i-2*step,0])
            ip1 = min([i+step,Np-1])
            ip2 = min([i+2*step,Np-1])
          END
          if keyword_set(linear) then $
                WaveTrans[i,j+1] = 0.5 * WaveTrans[i,j] + $
                              1./4. *  ( WaveTrans[im1,j]+WaveTrans[ip1,j] ) $
          else WaveTrans[i,j+1] = 3./8. * WaveTrans(i,j) + $
                           1./4. *  ( WaveTrans[im1,j]+WaveTrans[ip1,j] ) + $
                           1./16. * ( WaveTrans[im2,j]+WaveTrans[ip2,j] )
         i = i +1
      END
   WaveTrans[*, j] = WaveTrans[*, j] - WaveTrans[*, j+1]
END

DONE:

END
