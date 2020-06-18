function conversionfactor, FREQUENCIES,$
                           brightness2mKthermo=brightness2mKthermo,$
                           brightness2mKantenna=brightness2mKantenna,$
                           mKthermo2brightness=mKthermo2brightness,$
                           mKantenna2brightness=mKantenna2brightness,$
                           brightness2muKthermo=brightness2muKthermo,$
                           brightness2muKantenna=brightness2muKantenna,$
                           muKthermo2brightness=muKthermo2brightness,$
                           muKantenna2brightness=muKantenna2brightness,$
                           antenna2thermo=antenna2thermo,$
                           thermo2antenna=thermo2antenna
                           
;+
; NAME:
;	conversionfactor
;
; PURPOSE:
;       Delivers conversion factor between milli/mu K(antenna),
;       milli/mu K(thermodynamic) and brightness MJy/sr. Grand-unified subroutine.      
;
; CALLING SEQUENCE:
; 	A=conversionfactor(frequencies,/antenna2thermo)
;
; INPUTS:
;      FREQUENCIES REAL ARRAY    List of frequencies (GHz).
;
; OPTIONAL INPUTS:
;
; OPTIONAL INPUT KEYWORDS:
;      /brightness2mKthermo           Perform conversion MJy/sr -> mKthermo
;      /brightness2mKantenna          MJy/sr -> mKantenna
;      /mKthermo2brightness           mKthermo -> MJy/sr
;      /mKantenna2brightness          mKantenna -> MJy/sr
;      /brightness2muKthermo          Perform conversion MJy/sr -> muKthermo
;      /brightness2muKantenna         MJy/sr -> muKantenna
;      /muKthermo2brightness          muKthermo -> MJy/sr
;      /muKantenna2brightness         muKantenna -> MJy/sr
;      /antenna2thermo                Antenna -> thermodynamic
;      /thermo2antenna                Thermodynamic -> antenna
;
; NOTES
;
; SIDE EFFECTS
;
; EXAMPLES
;       Convert from antenna to thermodynamic units.
;
;       correctionfactor=conversionfactor(70,/antenna2thermo)
;
; COMMONS USED : 
;
; PROCEDURES USED: 
;
; MODIFICATION HISTORY:
; 	May 2006, Samuel Leach, SISSA
;-


;############Fundamental constants##############
c=2.99792458e10                 ; cm / s
h=6.626176e-27                  ; erg * s
k=1.380662e-16                  ; erg / K
T_cmb=2.726			; K

x= h*frequencies*1d9/k/T_cmb

;Just compute ALL conversion factors and be done with it.

conversion_thermo2antenna=x^2*exp(x)/(exp(x)-1.)^2
conversion_antenna2thermo=1./conversion_thermo2antenna

conversion_mKantenna2brightness=1.e-3*2.*k*(frequencies*1.e9/c)^2*1.e17
conversion_brightness2mKantenna=1./conversion_mKantenna2brightness
conversion_muKantenna2brightness=1.e-6*2.*k*(frequencies*1.e9/c)^2*1.e17
conversion_brightness2muKantenna=1./conversion_muKantenna2brightness

conversion_brightness2mKthermo=conversion_brightness2mKantenna*$
		   conversion_antenna2thermo
conversion_mKthermo2brightness=1./conversion_brightness2mKthermo
conversion_brightness2muKthermo=conversion_brightness2muKantenna*$
		   conversion_antenna2thermo
conversion_muKthermo2brightness=1./conversion_brightness2muKthermo


;Default output
conversion = conversion_antenna2thermo

if(keyword_set(brightness2mKthermo)) then conversion = conversion_brightness2mKthermo
if(keyword_set(brightness2mKantenna)) then conversion = conversion_brightness2mKantenna
if(keyword_set(brightness2muKthermo)) then conversion = conversion_brightness2muKthermo
if(keyword_set(brightness2muKantenna)) then conversion = conversion_brightness2muKantenna

if(keyword_set(mKthermo2brightness)) then conversion = conversion_mKthermo2brightness
if(keyword_set(mKantenna2brightness)) then conversion = conversion_mKantenna2brightness
if(keyword_set(muKthermo2brightness)) then conversion = conversion_muKthermo2brightness
if(keyword_set(muKantenna2brightness)) then conversion = conversion_muKantenna2brightness

if(keyword_set(antenna2thermo)) then conversion = conversion_antenna2thermo
if(keyword_set(thermo2antenna)) then conversion = conversion_thermo2antenna


return, conversion

END
