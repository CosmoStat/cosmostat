;+
; NAME:
;        survival
;
; PURPOSE:
;	 Return the survival function relative to a given distribution:
;         if nu < 0,     S(nu) =  #nbr_de_coeff / x < nu  
;         if nu >= 0,    S(nu) =  #nbr_de_coeff / x >= nu
;       If the keyword norm is set, the data are first normalized by:
;             Data = (Data- mean(Data)) / sigma(data)
;       If the keyword SigmaNorm is set, the data are first normalized by:
;             Data = Data / SigmaNorm
;       The nu threshold varies between -nsig and +nsig, and nsig has the value 10
;       by default. The keyword Np fixes the number of threshold values.
;       The output is a 2D IDL array:
;           Out[*,0]  = Nu values
;           Out[*,1]  = Survival values
;
; CALLING:
;
;      TabStat = survival( Data, norm=norm, SigmaNorm=SigmaNorm, Nsig=Nsig, Np=Np, plot=plot ) 
;       
; INPUTS:
;     Data -- IDL 2D array: Input data to analyze
;
; INPUT KEYWORDS:
;      Norm  : scalar -- if set, the input data are centered (i.e. Data = (Data-Mean)/Sigma)
;      SigmaNorm: scalar -- if set, the input data are normalized ((i.e. Data = Data/ SigmaNorm)
;      Nsig: float -- Maximum threshold value.
;      Np: int -- Number of thresholds
;      plot : scalar --  if set, plot  log( Out[*,1] + 1) versus Out[*,0]
;
; OUTPUT KEYWORDS: 
;
; EXAMPLE:
;       TabStat = survival(Data, /norm)
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	October, 2005 File creation
;-

;=========================================================================

function survival, Data, norm=norm, SigmaNorm=SigmaNorm, Np=Np, plot=plot, Nsig=Nsig
TabStat=0
if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE:  TabStat = survival(Data)'
        goto, DONE
        end
 	
N = N_ELEMENTS(Data)
X = dblarr(N)

if keyword_set(norm) then X[*] = (Data - mean(Data)) / sigma(data) $
else if keyword_set(SigmaNorm)then  X[*] = Data / SigmaNorm $
else X = Data

if not keyword_set(Nsig) then Nsig = 10.
if not keyword_set(Np) then Np = 1001
Np2 = Np / 2.

TabStat = dblarr(Np,2)
TabThreshold = 2. * dindgen(Np) / float(Np-1) * Nsig - Nsig  

TabStat [*,0] = TabThreshold 
for i=0,Np-1 do begin
   T = TabThreshold[i]
   if T lt 0 then begin 
      Ns = where(X LT T, c)
      if c  GT 0 then TabStat[i,1]=c else  TabStat[i,1]=0
   end else begin 
      Ns = where(X GT T, c)
      if c  GT 0 then TabStat[i,1]=c else  TabStat[i,1]=0
   end
end

if keyword_set(plot) then begin
      plot, TabStat[*,0], TabStat[*,1]  ;  ßalog(TabStat[*,1] + 1)
end
 

DONE:

return, TabStat

end

;=========================================================================


