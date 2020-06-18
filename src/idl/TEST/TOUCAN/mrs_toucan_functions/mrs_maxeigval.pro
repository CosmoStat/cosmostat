;+
; NAME:
;        mrs_maxeigval
;
; PURPOSE:
;       Estimate the maximum eigenvalue of a matrix by the power
;       iteration method.
;
; CALLING:
;     normT = mrs_maxeigval(Mat)
;
; INPUTS:
;     Mat        -- IDL array(M,M) = square matrix
;    
; OUTPUTS:
;     normT      -- IDL scalar = estimate of the maximum eigenvalue of
;                   the matrix Mat.
;
; EXAMPLE:
;      normT = mrs_maxeigval(Mat)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------

function mrs_maxeigval, Mat

Matsize = size(Mat)
M = Matsize[1]
normT = -1
if M ne Matsize[2] then begin
   print, 'Input matrix must be squared'
   goto, STOP
endif
tol = 0.01
difNorm = 1
x = randomn(seed, M)
x = x/sqrt(total(x^2))
temp = fltarr(M)
normT = 0
MaxIter = 5000
iter = 0
while ((difNorm gt tol) && (iter lt MaxIter)) do begin
   temp = Mat ## x
   difNorm = normT
   normT = sqrt(total(temp^2))
   difNorm = abs(difNorm - normT)
   x = temp/normT
   iter = iter + 1
   ;print, normT
endwhile
if iter eq MaxIter then print, 'max iteration reached'

STOP:
return, normT

end
