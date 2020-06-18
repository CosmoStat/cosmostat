;+
; NAME: 
;       DEL_PATTERN
;
; PURPOSE: 
;       Suppress in a pattern in an image.
;
; CALLING:
;       output = del_pattern(Imag, n_iter=n_iter, pattern=pattern,
;                            disp=disp, NbrScale=NbrScale,NSigma=NSigma)
;
; INPUTS:
;       Imag -- 2D IDL array: image in which we want to suppress a pattern
;
; KEYWORDS:
;       n_iter -- scalar: number of iterations (default is 3).
;       disp -- scalar: if set, then results are displayed during the iterations
;       NSigma -- float: Variance (default is 3.0)
;       NbrScale -- int: Number of scales used in the WT
;
; OUTPUTS:
;       pattern -- 2D IDL array: pattern = output - Imag
;
; HISTORY:
;	Written: Eric Pantin 1995.
;                Jean-Luc Starck: call iteratively remove_pattern - 1996
;-



;------------------------------------------------------------
function origin, frame
  ;; Defines the origin of even- and odd-sized frames.
  n = size(frame)
  result = intarr(2)
  result(0) = n(1)/2
  if (n(0) eq 2) then result(1) = n(2)/2 else result(1) = 0
  return, result
end
;------------------------------------------------------------
function centre, image
  ;; Shift the origin of the image from the bottom left corner to the centre.
  o = origin(image)
  s = size(image)
  if (s(0) eq 1) then $
    return, shift(image, o(0)) $
  else $
    return, shift(image, o(0), o(1))
end
;------------------------------------------------------------
function decentre, image
  ;; Shift the origin of the image from the centre to the bottom left corner.
  o = origin(image)
  s = size(image)
  if (s(0) eq 1) then $
    return, shift(image, -o(0)) $
  else $
    return, shift(image, -o(0), -o(1))
end
;------------------------------------------------------------
function dft, data
  ;; Discrete Fourier transform for centred origin pixel
  ;; The IDL FFT normalizes during the forward transformation.
  ;; return, fft(data, -1)
  return, centre(fft(decentre(data), -1))
end
;------------------------------------------------------------
function dfti, ft
  ;; Inverse discrete Fourier transform for centred origin pixel
  ;; return, fft(ft, +1)
  return, centre(fft(decentre(ft), +1))
end
;------------------------------------------------------------
pro med_mult_dec,image,I,J,nniveaux
dima=(size(image))(1)
dimb=(size(image))(2)
I=fltarr(dima,dimb,nniveaux+1)
J=fltarr(dima,dimb,nniveaux+1)

imalis_old=image
I(*,*,0)=image
for ii=1,nniveaux do begin

imalis_new=median(imalis_old,2*ii+1)
I(*,*,ii)=imalis_new
J(*,*,ii)=imalis_old-imalis_new
imalis_old=imalis_new
endfor
end

;------------------------------------------------------------

function remove_pattern, image_input, n_iter, disp=disp, NSigma=NSigma, NbrScale=NbrScale

image_work=image_input
image_ori=image_input
dima=(size(image_input))(1)
dimb=(size(image_input))(2)

if keyword_set(disp) then tvscl,congrid(image_input,256,256),0,256
if keyword_set(disp) then xyouts,40,260,'image entree',/device
if keyword_set(NbrScale) then Nb = NbrScale $
else Nb = 4

;estimation du sigma du bruit

;med_mult_dec,image_input,I,J,3

MROPT = '-n '+ STRCOMPRESS(STRING(Nb), /REMOVE_ALL)  

mr_transform, image_input, i, opt=MROPT
lastscale = Nb-1

;;on commence par enlever les structures a basse frequence
;; par median multiresolution

image_hf=image_ori-I(*,*,lastscale)
image_work=image_hf

sigma=stdev(image_hf,me_noise)

;;3 sigma clipping
for i=1,5 do begin
 ind_noise=where(abs(image_hf-me_noise) le 3.0*sigma,count_noise)
 ;verif=bytarr(dima,dimb)
 ;verif(ind_noise)=1
 ;display,verif
  if (count_noise ge 1) then sigma=stdev(image_hf(ind_noise),me_noise)  $
  else sigma=0
 
if keyword_set(disp) then  print,'new_sigma, moyenne',sigma,me_noise
endfor


;; remplacement des "objets" par une valeur mediane
;; un objet est detecte si il est au dessus de 3sigmas du bruit

ind_det_obj=where(abs(image_hf-me_noise) ge (3.0*sigma), count_obj)
if (count_obj ge 1) then image_work(ind_det_obj)=me_noise

for i_boucle=1,n_iter do begin
if keyword_set(disp) then  print,'boucle ',i_boucle
 image_work_ft=dft(image_work)
 image_work_ft_re=float(image_work_ft)
 image_work_ft_im=imaginary(image_work_ft)
if keyword_set(disp) then  tvscl,congrid(image_work_ft_re,256,256)
if keyword_set(disp) then  xyouts,40,20,'partie reelle fft',/device
if keyword_set(disp) then  tvscl,congrid(image_work_ft_im,256,256),256,0
if keyword_set(disp) then  xyouts,256+40,20,'partie imaginaire fft',/device

 ;;3 sigma clipping
 sigma_fl=stdev(image_work_ft_re,me_noise_fl)
 sigma_im=stdev(image_work_ft_im,me_noise_im)
 for i=1,5 do begin
  ind_noise_fl=where(abs(image_work_ft_re) le 3.*sigma_fl, count_fl)
  ind_noise_im=where(abs(image_work_ft_im) le 3.*sigma_im, count_im)
  if( count_fl ge 0) then sigma_fl=stdev(image_work_ft_re(ind_noise_fl),me_noise_fl)
  if( count_im ge 0) then sigma_im=stdev(image_work_ft_im(ind_noise_im),me_noise_im)
if keyword_set(disp) then   print,'sigmas float et imag',sigma_fl,sigma_im
 endfor

; image_work_ft_re( where(abs(image_work_ft_re) le $ 
;  (11-i_boucle)*sigma_fl) )=0.0
; image_work_ft_im( where(abs(image_work_ft_im) le $
;  (11-i_boucle)*sigma_im) )=0.0

image_work_ft_re( where(abs(image_work_ft_re) le NSigma*sigma_fl) )=0.0
image_work_ft_im( where(abs(image_work_ft_im) le NSigma*sigma_im) )=0.0

if keyword_set(disp) then  BEGIN
   tvscl,congrid(image_work_ft_re,256,256)
   xyouts,40,20,'partie reelle fft',/device
   tvscl,congrid(image_work_ft_im,256,256),256,0
   xyouts,256+40,20,'partie imaginaire fft',/device
END

 noise_pattern= float(dfti(complex(image_work_ft_re,image_work_ft_im)))
if keyword_set(disp) then  BEGIN
  tvscl,congrid(noise_pattern,256,256),512,0
  xyouts, 540,20,'noise pattern',/device
  tvscl,congrid(image_input-noise_pattern,256,256),256,256
  xyouts,280,280,'image nettoyee',/device
END
;; on masque les objets enleves et on les remplace par les pattern

if( count_obj ge 1) then image_work(ind_det_obj)=noise_pattern(ind_det_obj)
 
endfor


return,image_input-noise_pattern

end
;------------------------------------------------------------

function del_pattern, image_input, n_iter=n_iter, pattern=pattern, disp=disp, NSigma=NSigma, NbrScale=NbrScale

if N_PARAMS() LT 1 then begin 
        print, 'CALL SEQUENCE: output = del_pattern(Imag, n_iter=n_iter, pattern=pattern, disp=disp, NSigma=NSigma)'
        goto, DONE
        end

if keyword_set(disp) then window,0,xsize=256*3,ysize=512
if not keyword_set(NSigma) then NSigma = 3.
if not keyword_set(NbrScale) then NbrScale = 4
print, "Number of scales = ", NbrScale

if not keyword_set(n_iter) then n_iter = 3
data = image_input
for i=0, n_iter do BEGIN
    image_output = remove_pattern(data,1,disp=disp,NSigma=NSigma,NbrScale=NbrScale)
    pattern = image_input - image_output
END
return, image_output

DONE:

END

