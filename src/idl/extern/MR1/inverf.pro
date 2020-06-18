;+
; NAME:
;INVERF -- inverse error function
;     
; PURPOSE:
;       Returns the inverse error function for a given probabality.
;       
;     
; CALLING SEQUENCE:
;       RESULT = INVERF(X)
;     
; INPUTS:
;       X : Number representing the integrated Gaussian probability 
;           distribution between -N standard deviations and +N standard
;           deviations.
;     
; OUTPUTS:
;       Returns the inverse error function of X.  This is the number of
;       standard deviations that the error function needs to be integrated 
;       over in order to yield the value X.
;
; KEYWORDS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       The error function is normalized to unity, so it is meaningless
;       to ask for the inverse error function of a value larger than
;       one or less than zero.  If you do, the function will stop and 
;       tell you to contemplate your mistake.
;
;       To check behavior visually:
;         IDL> x = dindgen(10001)/10000.
;         IDL> plot, x, abs(x-errorf(inverf(x))), /ylog, yr=[1d-9,1d-1]
;         IDL> oplot, !x.crange, 4.5d-4*[1,1], lines=1
;
;       Hasting's approximation (see NOTES below) is accurate to
;       4.5d-4, and we see this is indeed the case.
;
;       Alternatively, check how large x can become before
;       inverf(errorf(x)) fails, i.e. the difference between
;       x and inverf(errorf(x)) becomes larger than 4.5d-4:
;
;         IDL> x = dindgen(10001)/10000.*6
;         IDL> plot, x, abs(x-inverf(errorf(x))), /ylog, yr=[1d-9,1d-1]
;         IDL> oplot, !x.crange, 4.5d-4*[1,1], lines=1
;
;       So, if x is double precision, inverf(errorf(x)) is accurate
;       up to x=5.34.  If x is floating precision, it is accurate up to
;       x=3.11.  This makes plenty of sense, since beyond this limit out in 
;       the wings of the normal distribution, the error function is exactly 
;       one to machine precision. So the inverse error function can't be
;       expected to differentiate between the error function of 7 and the
;       error function of 200!
;
; NOTES:
;       Using notation from Abramowitz & Stegun, page 931, we see:
;
;         erf(x) = 2*P(x*sqrt(2)) - 1
;         P(x*sqrt(2)) = 1 - Q(x*sqrt(2))
;
;       Given Q(x*sqrt(2))=p, we can use Hasting's approximation 
;       for digital computers (Abramowitz & Stegun, page 933) to get:
;
;         x*sqrt(2) = t - (c0+c1*t+c2*t^2)/(1+d1*t+d2*t^2+d3*t^3) + e(p)
;
;       where |e(p)| < 4.5d-4 and the constants appear in the code below.
;
; MODIFICATION HISTORY:
;       Written Tim Robishaw, Berkeley  22 Feb 2002
;       Stole idea from Carl Heiles.
;       Added machine precision considerations. TR 23 Feb 2002
;-

function inverf, erfx

on_error, 2

; ERROR FUNCTION CAN ONLY RETURN 0 <= ERRORF(X) <= 1, SO
; THE INVERSE ERROR FUNCTION OF ANYTHING GREATER THAN ONE IS
; MEANINGLESS...
if (total((erfx le 1d0) AND (erfx ge 0)) ne N_elements(erfx)) then $
    message, 'The range of the error function is 0 < erf(x) <= 1!'

; IS ERFX WITHIN MACHINE PRECISION OF UNITY...
type = size(erfx,/TYPE)
epsn = (machar(DOUBLE=(type eq 5L))).epsneg

; CHECK THE MACHINE PRECISION HERE...
p = 1d0 - (erfx - epsn*( (erfx - 0.5*epsn) eq fix(1.,TYPE=type) ))
p = 0.5d0*p
t = sqrt(-2d0*alog(p))
num =  (0.010328d0*t + 0.802853d0)*t + 2.515517d0
den = ((0.001308d0*t + 0.189269d0)*t + 1.432788d0)*t + 1d0
return, 0.70710678118654752440d0 * ((t - num/den) > 0)

end; inverf
