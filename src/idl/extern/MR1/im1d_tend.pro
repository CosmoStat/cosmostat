;+
; NAME:
;        IM1D_TEND
;
; PURPOSE:
;	Remove the tendancy (baseline) in a 1D signal using a polynomial fit.
;
; CALLING:
;
;      IM1D_TEND, SigIn, SigOut, Tend=Tend, WinSize=WinSize
;
; INPUTS:
;     Sig -- 1D IDL array: Input Signal
;
; OUTPUTS:
;     SigOut -- 1D IDL array: Output signal
;
; KEYWORDS:
;      Tend -- 1D IDL array (outpout): estimated base line
;
;      WinSize --  integer: window size for the polynomial fit.
;                           Default is 100.
;
; EXTERNAL CALLS:
;       im1d_tend (C++ program)
;
; EXAMPLE:
;
;       Remove the baseline using a window size of 60 pixels
;              RES = IM1D_TEND(Sig1d, WinSize=60)
;
; HISTORY:
;	Written: Jean-Luc Starck 2004.
;	July, 2004 File creation
;-

pro IM1D_TEND,  SigIn, SigOut, Tend=Tend, WinSize=WinSize

SigIn = float(SigIn)
if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: IM1D_TEND, SigIn, SigOut, Tend=Tend, WinSize=WinSize'
        goto, DONE
        end

vsize = size(SigIn)
if vsize(0) NE 1 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM1D_TEND, SigIn, SigOut, Tend=Tend, WinSize=WinSize'
         goto, DONE
        end
 
if not keyword_set(Opt) then Opt = ' '

filename = 'xx_temp.fits'
 p = strpos(filename, '.fits')
if p LT 0 then filename = filename + '.fits'

NameSig='xx_imag.fits'
NameTend='xx_tend.fits'
writefits, NameSig, SigIn

if  keyword_set(WinSize) then  opt= ' -T ' + STRCOMPRESS(STRING(WinSize), /REMOVE_ALL) $
else opt=' '

com = 'im1d_tend ' + ' ' + Opt + ' ' + NameSig + ' ' + NameTend + ' ' +  filename
spawn, com
SigOut = readfits(filename)
Tend = readfits(NameTend)

delete, NameSig
delete, NameTend
if filename EQ 'xx_temp.fits' then delete, filename
DONE:

END
