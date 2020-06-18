pro get_prt_temp, serial_no, resist, temperature, Sim_PRT=sim_prt
;+
; NAME:
;      GET_PRT_TEMP
;
; PURPOSE:
;      Convert PRT (platinum resistance thermometer) resistances in ohms to 
;      temperature in Kelvin. User must know the serial number of the PRT of 
;      interest.
;
; CALLING SEQUENCE:
;      get_prt_temp, serial_no, resist, temperature
;
; INPUTS:
;      serial_no - string - serial number of the PRT, e.g., 'UG84'.
;                           Case insensitive.
;                           Note that only ONE serial number is
;                           allowed per calling sequence!
;      resist    - fltarr - an array of measured resistances for
;                           the specified PRT, in ohms. The input array
;                           may contain one data point or many data points.
;                           Resistance values are expected to lie within
;                           the range of approximately 10 - 650 ohms.
;
; OUTPUTS:
;      temperature - fltarr - temperatures, in Kelvin, which correspond
;                             to the input resistances.
;                             Temperature is set to -1 if an input
;                             resistance lies outside the valid interpolation
;                             range.
;
; OPTIONAL INPUT KEYWORDS:
;      /Sim_PRT - If present and nonzero, used the fake conversion curve
;                from fake_ihk to convert from resistance to K.
;
; COMMENTS:
;      Algorithm used here is based on the fortran routine find_t.f90
;      from S. Meyer/J. Grimes.  The IDL procedure returns identical
;      results to that of the fortran.
;      Temperatures are interpolated using b-splines.  B-spline coefficients
;      are restored from an IDL saveset, prt_splinedata.xdr, located in 
;      the reference area.  Data contained in this saveset originated from
;      an ASCII file provided by S. Meyer/J. Grimes.
;     
; EXAMPLE:
;      IDL> get_prt_temp,'UG84',[300.,400.],temps
;           ===> temps = [174.171, 223.509] 
;
; MODIFICATION HISTORY:
;      Original version: J. Weiland,  07 January 1999
;      Sim_PRT keyword added.  JW, June 2000.
;      Use getenv to extract environment variables.  MRG, SSAI, 20 August 2002.
;-
;-------------------------------------------------------------------------
;
; Conversion for fake ihk prt data- up top
;
IF keyword_set(sim_prt) then begin
  temperature = (resist + 68.d0)/2.06154d0
  return
endif

; For REAL data...

;
; Retrieve the spline data coefficients for each PRT
;
restore,concat_dir( getenv('MAP_REF'), $
  concat_dir('prt_data','prt_splinedata.xdr') )
;
; Find the correct data structure element for this serial no
;
serial_no = strtrim(strupcase(serial_no),2)
element = where (strtrim(string(splinedata.serial),2) eq serial_no, Nfound)
if Nfound GT 0 then begin
    element = element[0]
endif else begin
    print,'Serial number ',serial_no,' not recognized! Recheck and try again.'
    return
endelse
;
; Get the spline structure data for this PRT only
;
a = splinedata[element]
;
; Giant loop for the resistance arrary
;
n_resist = n_elements(resist)
temperature = fltarr(n_resist)

for jp = 0L, n_resist-1  do begin

  resistance = resist[jp]
;
; Find the index j such that xknot(j) .le. resistance .le. xknot(j+1)
;
  indices = where ((a.xknot lt resistance) and (a.xknot ne 0.0)) 
  j = max(indices)
;
; Protect against extrapolation
;
  minknot = min(a.xknot[where(a.xknot ne 0.0)])
  maxknot = max(a.xknot)

  if ((resistance lt minknot) or (resistance gt maxknot)) then begin
     temperature[jp] = -1.0

  endif else begin                   ; interpolation OK to proceed
;
; initialize b-spline coefficients
;
     k = a.korder
     d = fltarr(24,8)                            ; maximum workspace necessary
     d[0:a.ncoef-1,0] = a.bscoef[0:a.ncoef-1]
;
; find value of b-spline using de Boor algorithm
;
     for r = 1, k-1 do begin
       for i = (j-(k-1)),j-r  do begin
           d[i,r] = ((a.xknot[i+k] - resistance)*d[i,r-1])     $
                      / (a.xknot[i+k] - a.xknot[i+r])          $
                    + ((resistance-a.xknot[i+r])*d[i+1,r-1])   $
                      /(a.xknot[i+k] - a.xknot[i+r])
       endfor
     endfor        

     temperature[jp] = d[j-(k-1),k-1]

  endelse
endfor
;
;
return
end
