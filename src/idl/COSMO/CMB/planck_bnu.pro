FUNCTION PLANCK_BNU, nu, temp=temp
;+
; NAME:
;      PLANCK_BNU
; PURPOSE:
;      Compute Planck blackbody spectrum (W/m2/sr/Hz) at given
;      frequencies nu for a given temperature temp (CMB temperature 
;      of !const.tcmb by default)    
; CALLING SEQUENCE:
;      result = PLANCK_BNU(nu, temp=temp)
; INPUT:
;      nu : frequencies in Hertz, may be an array
; OUTPUT:
;      Planck spectrum at frequencies nu
; KEYWORDS:
;      temp: temperature of the blackbody
; RESTRICTIONS:
;       system variable giving physical constant !const.tcmb needed
;       Not fully tested yet. Report any bugs to Jacques Delabrouille
;       j.delabrouille@cdf.in2p3.fr
; PROCEDURES CALLED:
;       none
; REVISION HISTORY
;       Written, Jacques Delabrouille July 1999
;	Name changed from JD_BNU to PLANCK_BNU by JD, February 2006
;-

c = 299792458d
h = 6.62607554d-34
k = 1.38065812d-23
IF KEYWORD_SET(temp) EQ 0 THEN temp = !const.tcmb 
	;= 2.725 K... the system variable !const should have been defined
	;prior to using this routine
RETURN, 2*h*nu*(nu/c)^2 / (EXP(h*nu/k/temp) - 1d)
END