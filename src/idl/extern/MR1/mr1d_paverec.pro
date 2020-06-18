;+
; NAME:
;        MR1D_PAVEREC
;
; PURPOSE:
;	Reconstruct a signal from its pyramidal multiresolution transform
;
; CALLING:
;
;      MR1D_PAVEREC, MR_Trans, Rec_Signal, Adjoint=Adjoint, 
;                     mirror=mirror, nosmooth=nosmooth
;
; INPUTS:
;     MR_Trans -- 2D IDL array:  multiresolution transform without undersampling
;
; OUTPUTS:
;     Rec_Signal -- 1D array: reconstructed signal
;
; KEYWORDS
;     Adjoint -- scalar: if set, the adjoint operator is used for the 
;                        reconstruction
;     nosmooth -- scalar: if set, the last scale is not used.
;
;     mirror -- scalar: if set, and if Adjoint is set, 
;                        then a mirror is used at the border
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-

pro smooth1d, Signal, DataSmooth, step=step, mirror=mirror, linear=linear

Np = (size(Signal))[1]
DataSmooth = Signal
DataSmooth[*]=0
if not keyword_set(step) then step=1
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
              DataSmooth[i] = 0.5*Signal[i] + 1./4.*(Signal[im1]+Signal[ip1]) $
          else DataSmooth[i] = 3./8. * Signal[i] + $
                           1./4. *  ( Signal[im1]+Signal[ip1] ) + $
                           1./16.*  ( Signal[im2]+Signal[ip2] )
         i = i +1
   END
End

;-----------------------------------------------------

pro MR1D_PAVEREC, MR_Trans, Rec_Signal, Adjoint=Adjoint, mirror=mirror, nosmooth=nosmooth

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: MR1D_PAVEREC, MR_Trans, Rec_Signal, Adjoint=Adjoint, mirror=mirror'
        goto, DONE
        end

Nscale = (size(MR_Trans))[2]
Rec_Signal = MR_Trans[*, Nscale-1]
if keyword_set(nosmooth) then Rec_Signal[*] = 0

for j=Nscale-2,0,-1 do BEGIN
    if keyword_set(Adjoint) then BEGIN
      if j GT 0 then BEGIN
             step = 2^(j-1)
             smooth1d, Rec_Signal, Smooth, step=step, mirror=mirror
             Rec_Signal = Smooth
;             step = 2^j
;             smooth1d, MR_Trans[*, j], Smooth, step=step, mirror=mirror
;             Rec_Signal = Rec_Signal + MR_Trans[*, j] - Smooth
             Rec_Signal = Rec_Signal + MR_Trans[*, j]
;             if j NE Nscale-2 then BEGIN
;                step = 2^j
;                smooth1d, MR_Trans[*, j], Smooth, step=step, mirror=mirror
;               Rec_Signal = Rec_Signal + MR_Trans[*, j] - Smooth
;             END ELSE Rec_Signal = Rec_Signal + MR_Trans[*, j]
      END ELSE BEGIN
              step = 2^j
             smooth1d, MR_Trans[*, j], Smooth, step=step, mirror=mirror
             Rec_Signal = Rec_Signal + MR_Trans[*, j] - Smooth
;             Rec_Signal = Rec_Signal + MR_Trans[*, j]
      END
    END else Rec_Signal = Rec_Signal + MR_Trans[*, j]
   END


DONE:

END
