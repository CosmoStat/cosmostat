;+
; NAME:
;        MR1D_DETECT_LEVEL
;
; PURPOSE:
;	Separate high and low frequencies, and gives the noise level
;       for each band.
;
; CALLING:
;
;      MR1D_DETECT_LEVEL, Signal, Nscale, Wt_Signal, Wt_Level, 
;                         continuum=continuum, median=median,
;                         LFreq_Level=LFreq_Level, H_Freq=H_Freq, L_Freq=L_Freq,
;                         Error_H_Freq=Error_H_Freq, Error_L_Freq=Error_L_Freq
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want the continuum
;    
;     Nscale -- scalar: number of scales
;
; INPUT KEYWORDS:
;
;     median=median -- scalar: if set, use the pyramidal median transform,
;                               (and not the a trou wavelet transform).
;
;     LFreq_Level -- scalar: high frequencies are from 0 to LFreq_Level-1
;                            and low frequencies are from LFreq_Level to Nscale
;                            This parameter is used for H_Freq and L_Freq
;                            estimation.
;
;     mirror -- scalar: if set, and if a trous wavelet transform is used,
;                       the borders are manadged by mirror effect.
;
; OUTPUTS:
;     Wt_Signal -- 1D IDL array: Wavelet transform of the signal
;
;     Wt_Level -- 1D IDL array: Noise estimation at each scale
;
;
; INPUT-OUTPUT KEYWORDS:
;
;     continuum -- 1D IDL array: continuum of the spectrum
;     TabSigmaNoise -- 1D IDL array: RMS of the input data
;
; OUTPUT KEYWORDS:
;
;     H_Freq -- 1D IDL array: high frequencies (from 0 to LFreq_Level-1)
;     L_Freq -- 1D IDL array: low frequencies  (from LFreq_Level to Nscale)
;     Error_H_Freq -- 1D IDL array: noise upper limit on high frequencies
;     Error_L_Freq -- 1D IDL array: noise upper limit on low frequencies
;     
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR1D_DETECT_LEVEL, Signal, Nscale, Wt_Signal, Wt_Level,  $
                     continuum=continuum, median=median, $
                     LFreq_Level=LFreq_Level, H_Freq=H_Freq, L_Freq=L_Freq, $
                     Error_H_Freq=Error_H_Freq, Error_L_Freq=Error_L_Freq, $
                     TabSigmaNoise=TabSigmaNoise, mirror=mirror

tau=-1
if N_PARAMS() LT 2 then begin 
      print, 'CALLING SEQUENCE: MR1D_DETECT_LEVEL, Signal, Nscale, '
      print, '                    Wt_Signal, Wt_Level,continuum=continuum,'
      print, '                    median=median, LFreq_Level=LFreq_Level,'
      print, '                    H_Freq=H_Freq, L_Freq=L_Freq,'
      print, '                    Error_H_Freq=Error_H_Freq,'
      print, '                    Error_L_Freq=Error_L_Freq'
      goto, DONE
      end

Np = (size(Signal))(1)
if Np LT 2 then  begin 
      print, 'bad first parameter ...'
      goto, DONE
      end

if not keyword_set(continuum) then continuum = mr1d_continuum(Signal, Nscale)

 ; standard deviation of the noise at each sample 
NoiseSignal = Signal - continuum
if not keyword_set(TabSigmaNoise) then BEGIN
   TabSigmaNoise = fltarr(Np)
   for k=0, Np-1 do $
   BEGIN
     if k lt 10 then TabSigmaNoise(k) = sigma(NoiseSignal(0:20)) $
     else BEGIN
        if k GT Np-11 then TabSigmaNoise(k) = sigma(NoiseSignal(Np-11:Np-1)) $
        else TabSigmaNoise(k) = sigma(NoiseSignal(k-10: k+10))
     END
   END
END

; detection level at each scale
TabLevel = [0.726336,0.283629,0.177784,0.1261,0.0846,0.0597747, 0.0386936, 0.0307165, 0.0223076]
TabLevelMed = [0.915152,0.357905,0.248778,0.209270,0.162813,0.111927, 0.0786654, 0.0594114, 0.0643011]
TabSmooth = [1.,0.521788,0.352123,0.244073,0.167545, 0.116704, 0.0816195, 0.0816195, 0.0816195]
TabSmoothMed = [1.,0.531919,0.399170,0.311026,0.239450, 0.178777, 0.139950, 0.139950, 0.139950]

if keyword_set(median) then BEGIN
   TabLevel=TabLevelMed
   TabSmooth=TabSmoothMed
  ; mr1d_pyrmed, NoiseSignal, Wt_Signal, Nscale, /interp
     mr1d_pavemed, NoiseSignal, Wt_Signal, Nscale
END ELSE  BEGIN
 if keyword_set(mirror) then mr1d_atrou,NoiseSignal,Wt_Signal,Nscale,/mirror $
 else mr1d_atrou,NoiseSignal,Wt_Signal,Nscale
END

Wt_Level=Wt_Signal
Wt_Level(*)=0
wt_err=Wt_Level
for j = 0, Nscale-1 do $
BEGIN
      WinSize = 2^(j+1)
      for k=0, Np-1 do BEGIN
         Ma = max([0,k-WinSize])
         Mi = min([Np-1,k+WinSize])
         wt_err(k,j) =  max(TabSigmaNoise(Ma:Mi))
         Wt_Level(k,j) =  wt_err(k,j)*TabLevel(j)
      END
END

; low and high frequencies estimation
if keyword_set(LFreq_Level) then BEGIN
   H_Freq = fltarr(Np)
   L_Freq = fltarr(Np)
   H_Freq_Error = fltarr(Np)
   L_Freq_Error = fltarr(Np)
   for f = 0,LFreq_Level-1 do H_Freq = H_Freq + Wt_Signal(*,f) 
   for f = LFreq_Level,Nscale-1 do L_Freq = L_Freq + Wt_Signal(*,f)
   Error_H_Freq = wt_err(*, LFreq_Level-1)*TabLevel(0)
   Error_L_Freq = wt_err(*, NScale-1)*TabSmooth(LFreq_Level)
END
DONE:
return
END
