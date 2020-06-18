;-------------------------------------------------------------------------------
;
; NAME:
;        ATWT1D
;
; PURPOSE:
;       Computes the a trous wavelet transform of a signal. 
;       The output is a 2D IDL array. 
;
; CALLING:
;
;      ATWT1D, Signal, TransSignal, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord
;       
; INPUTS:
;     Signal -- 1D IDL array: signal we want to transform
;
; OUTPUTS:
;     TransSignal -- 2D IDL array: Wavelet Transform
;                                 TransSignal(i, *) = ith band of the 
;                                                    wavelet transform
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
;                     TransSignal is an array of Nscale signals: 
;                                      TransSignal(0: Nscale-1, *)
;			There is no test of the validity of the number of scales with respect
;			to the size of the input images. 
;
;      TabStat -- 2D IDL Array: Statistic of the transform
;                               TabStat[i,0] = mean value of band i
;                               TabStat[i,1] = sigma value of band i
;                               TabStat[i,2] = skewness value of band i
;                               TabStat[i,3] = kurtosis value of band i
;
;      NoBord -- int: if NoBord is set, the wavelet coefficients 
;                     which present a border problem are not taken into
;                     account in the statistics calculation.
; 
; EXAMPLE:
;
;       Compute the multiresolution of the signal I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               ATWT1D, I, Output
;
; HISTORY:
;       Written: Jean-Luc Starck 2002.
;       June, 2002 File creation
;
; REFERENCE:
;    J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
;
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
function test_ind_1D, ind, N

	;this function implements mirror like limit conditions on the edges of the input signal
	;ATTENTION : the output may still be out of range ie not in [0, N-1]
	;refinements may be necessary, although not meaningful in practice
	
	ret = ind
	if ind LT 0 then ret = -ind $
	else if ind GE N then ret = 2*N-ind-2
	
	return, ret
	
end
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
pro b3spline_1D, Sig_in, Sig_out,  Step
   
   N = double( (size(Sig_in))(1) )
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Sig_out  = Sig_in

   for i = 0., N-1. do begin
   
      im = test_ind_1D(i-Step, N)
      ip = test_ind_1D(i+Step, N)
      im2 = test_ind_1D(i-2*Step, N)
      ip2 = test_ind_1D(i+2*Step, N)
      Sig_out(i) = C3 * Sig_in(i) + C2 * (Sig_in(im) + Sig_in(ip)) + C1 * (Sig_in(im2) + Sig_in(ip2)) 

   endfor

end
;-------------------------------------------------------------------------------




;-------------------------------------------------------------------------------
pro atwt1d, Signal, TransSignal, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord


	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: atwt1d, Signal, TransSignal, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4

	if Nscale LT 2 then  Nscale =4
	
	;-------------------------------------------------------------------------------
	;recovering a few parameter values, initializing the transform loop ...
	N = (size(Signal))(1)

	NStep=Nscale-1
	Nstat = 4
	TabStat= fltarr(Nstat,Nscale)
	TransSignal = fltarr(Nscale,N)
	Signal_in = Signal
	Signal_Out = 0.* Signal
	Step_trou = 1

	for i=0., NStep-1 do begin
		;B3 Spline smoothing
		b3spline_1d, Signal_in, Signal_Out,  Step_trou

		; Wavelet coefficient calculations
		TransSignal(i, *) = Signal_in - Signal_Out
		Signal_in = Signal_Out

		; Statistic calculation
		BordSize = 2*Step_trou
		if not keyword_set(NoBord) then BordSize = 0
	
		band = TransSignal(i, BordSize:N-1-BordSize )
		TabStat[0,i] = total(band, /double) / N_ELEMENTS(band)
		TabStat[1,i] = sigma(band)  
		TabStat[2,i] = skewness(band) 
		TabStat[3,i] = kurtosis(band) 
		;print, "Scale ", i+1, "  mean = ",  TabStat[0,i] ,"   Sigma = ", TabStat[1,i], "   Skewness = ", TabStat[2,i], "   Kurtosis = ", TabStat[3,i]
   
		; New distance between two pixels
		if i NE NStep-1 then Step_trou = Step_trou * 2.
   endfor


	; Smooth array
	i =  NStep  
	TransSignal(i,*) = Signal_Out
	band = TransSignal(i, BordSize:N-1-BordSize )
	TabStat[0,i] = total(band, /double) / N_ELEMENTS(band)
	TabStat[1,i] = sigma(band)  
	TabStat[2,i] = skewness(band) 
	TabStat[3,i] = kurtosis(band) 
	;print, "Smooth Ima:  ", "  mean = " , TabStat[0,i] , "    Sigma = ", TabStat[1,i], "   Skewness = ", TabStat[2,i], "   Kurtosis = ", TabStat[3,i]


DONE:

end


;-------------------------------------------------------------------------------

