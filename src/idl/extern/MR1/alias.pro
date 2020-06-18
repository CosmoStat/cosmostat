
;;;----------------------------------------------------------------------------
;;; Get out of IDL by all possible means.
;;;----------------------------------------------------------------------------

;------------------ QUIT --------------

pro mrh
widget_olh
end

pro mrsh
widget_olh
end

pro quit
  ;; An alias for IDL's EXIT command.
  exit
  return
end

pro logout
  ;; An alias for IDL's EXIT command.
  exit
  return
end

pro bye
  ;; An alias for IDL's EXIT command.
  exit
  return
end

; ---------------- delete -------------

pro delete, filename
case !version.os of
    'vms': cmd = 'delete'
    'windows': cmd = 'del'
    'MacOS': goto, notsup
    else: cmd = '\rm -f'
    endcase
cmd = cmd + ' ' + filename  
spawn, cmd  
goto, DONE
NOTSUP: print, "This operation is not supported on the Macintosh'
DONE:
end

;------------------ LS ----------------

pro ls
case !version.os of
    'vms': cmd = 'dir'
    'windows': cmd = 'dir'
    'MacOS': goto, notsup
    else: cmd = 'ls'
    endcase
if !version.os NE 'vms' then spawn, cmd, /noshell   $
else  spawn, cmd
goto, DONE
NOTSUP: print, "This operation is not supported on the Macintosh'
DONE:
end

;--------------- HISTORY----------------

pro h
help, /RECALL_COMMANDS
end

pro hs, x
help, /struct, x
end

pro hh, x
help, x
end 

;-------------- MEMORY------------------
pro pstat
help, /MEMORY
end
;---------------- PWD ------------------

pro pwd
case !version.os of
    'vms': cmd = 'show def'
    'windows': cmd = 'pwd'
    'MacOS': goto, notsup
    else: cmd = 'pwd'
    endcase
if !version.os NE 'vms' then spawn, cmd, /noshell $
else spawn, cmd
goto, DONE
NOTSUP: print, "This operation is not supported on the Macintosh'
DONE:
end

;--------------- LUT -----------------

pro lut, i
  ;; An alias for IDL's LOADCT command.
  loadct, i
  return
end

;pro wlut
;xloadct, /USE_CURRENT
;end

pro luts
  ;; Loop through all look-up tables. 
  ;; Actually this is not an alias but a new command.
  for i = 0, 15 do begin
    loadct, i
    wait, 2
  endfor
  return
end

;------------- min, max, sgn -----------------

function my_max, a, b
  ;; Return maximum of two numbers.
  if (a gt b) then begin
    return, a
  endif else begin
    return, b
  endelse
end

function my_min, a, b
  ;; Return minimum of two numbers.
  if (a lt b) then begin
    return, a
  endif else begin
    return, b
  endelse
end

function sgn, x
  ;; Returns the sign of an array.
  result = 0*x
  pos_supp = where(x gt 0, pos_cnt)
  if (pos_cnt gt 0) then result[pos_supp] = +1
  neg_supp = where(x lt 0, neg_cnt)
  if (neg_cnt gt 0) then result[neg_supp] = -1
  return, result
end


FUNCTION SIGN, X
return, fix(x gt 0.) - fix(x lt 0.)
end

;;;----------------------------------------------------------------------------
;;; Image display.
;;;----------------------------------------------------------------------------
pro tvilut, image, window=window, title=title, horizontal=horizontal, nowin=nowin
   if  not keyword_set(window) then w = 0 else w = window
   sz=size(image)
   if not keyword_set(nowin) then begin
      if not keyword_set(horizontal) then window, w , xsize=sz(1)*1.2, ysize=sz(2), title=title $
      else window, w , xsize=sz(1), ysize=sz(2)*1.2, title=title
   end else clear
   tvscl, image
   if not keyword_set(horizontal) then colorbar, /vertical, MINRANGE=min(image), MAXRANGE=max(image), format='(F18.1)',position= [.95, 0.10, .98, 0.90] $
   else colorbar,  MINRANGE=min(image), MAXRANGE=max(image), format='(F18.1)',position= [.95, 0.10, .98, 0.90]
end

pro display, image
  ;; An alias for IDL's TVSCL command.
  tvscl, image
  return
end

pro disp, window=window
if  not keyword_set(window) then w = 0 else w = window
window, w, xsize=512, ysize=512
  return
end


pro winb, win=win
if not keyword_set(win) then window,  xsize=1024, ysize=1024 $
else window, win, xsize=1024, ysize=1024
end

pro winbs, win=win
if not keyword_set(win) then win=0
window, win, xsize=1600, ysize=500
end

pro tvi, ima, win=win, set=set, size=size
if keyword_set(size) then begin
  ima1 = congrid(ima, size, size)
  vs=size(Ima1)
end else vs=size(Ima)
nx=vs[1]
ny=vs[2]
if keyword_set(win) then window, win, xsize=nx, ysize=ny
if keyword_set(set) then wset, set

if keyword_set(size) then tvscl, ima1 $
else tvscl, ima
end

pro tvis, ima, window=window 
vs=size(Ima)
nx=vs[1]
ny=vs[2]
if keyword_set(window) then window, window, xsize=nx, ysize=ny $
else window, xsize=nx, ysize=ny
tvscl, ima
end

function bgdisp, dat

n = size(dat) & Nc = n[1] & Nl = n[2]
Coefx = 512. / float(Nc)
Coefy = 512. / float(Nl)

Coef = my_min (Coefx, Coefy)
if Coef LT 1. then Coef =  my_max (Coefx, Coefy)
Nlb = Nl*Coef
Ncb = Nc*Coef
return, congrid(dat, Ncb, Nlb)
end

pro load, image 
tvscl, bgdisp(image)
 return
end

pro frame, cube, i
  a = cube[*,*,i]
  a = reform(a,/overwrite)
  tvscl, congrid(a, 256,256) 
end

pro diso, im
  tvscl, congrid(im, 256,256)
  end

pro load_cube, cube
   imax = (size(cube))[3]
   window, 1, xsize=256, ysize=256
   for i = 0, imax-1 do begin
    print, i 
    a = cube[*,*,i]
    a = reform(a,/overwrite)
    tvscl, congrid(a, 256,256)
;    load, a
  endfor
 end
 
pro clear
  ;; An alias for IDL's erase command.
  erase
  return
end

pro fblink, a, b, rate
  ;; An alias for IDL User Library's FLICK command.
  ; flick, a, b, rate
  flick, a*(255/max(a)), b*(255/max(b)), rate
  return
end

pro spec, s
plot, s, /YNOZERO
end

;;;----------------------------------------------------------------------------
;;; Evaluator
;;;----------------------------------------------------------------------------
pro eval, statement
  ;; Evaluates an IDL statement. An alias for IDL's EXECUTE function.
  e = execute(statement)
  return
end
;; eval, "help"

;;;----------------------------------------------------------------------------
;;; Miscellaneous.
;;;----------------------------------------------------------------------------

pro pause, seconds
  ;; An alias for IDL's WAIT command.
  wait, seconds
  return
end

pro abort
  ;; An alias for IDL's retall
  retall
  return
end

pro top
  ;; An alias for IDL's RETALL command.
  retall
  return
end

pro man, cmd
doc_library, cmd
end

;;;----------------------------------------------------------------------------
;;; Math functions
;;;----------------------------------------------------------------------------

function re, x
  return, float(x)
end

function im, x
  return, imaginary(x)
end

function mean, array
  ;; Return the average of a 1D or 2D data array.
  npix = size(array)
  return, total(array) /  N_ELEMENTS(array)
end

pro info, Tab
print, "Min = ", min(Tab)
print, "Max = ", max(Tab)
print, "mean = ", mean(Tab)
print, "sigma = ", sigma(Tab)
print, "total = ", total(Tab)
end


;==============================

function mygauss1d, idim, sigma, pos=pos

;-----------------
; parameters check
;-----------------
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, $
   'CALLING SEQUENCE: GAUSS=  mygauss1d( idim, sigma, pos=pos)'
   GOTO, CLOSING
 ENDIF
if  not keyword_set( pos) then begin
   pos = 0
 endif
x = findgen(idim) - idim/2.

GAUSS= EXP(-(x - POS )^2 /(2.*SIGMA^2))

return, gauss
closing:
end

;==============================

function mygauss, idim, jdim, sigma, pos=pos

;-----------------
; parameters check
;-----------------
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, $
   'CALLING SEQUENCE: GAUSS=  mygauss( idim, jdim, sigma, pos=pos)'
   GOTO, CLOSING
 ENDIF
if n_elements( pos) eq 0 then begin
  pos= make_array( 2, /float)
  pos(0)= FLOAT(IDIM)/2.0
  pos(1)= FLOAT(JDIM)/2.0
endif
if n_elements( pos) eq 1 then begin
  print, 'Give 2 numbers for POS'
  goto, closing
endif

GAUSS= EXP(-(( FINDGEN( IDIM)#REPLICATE( 1.0, IDIM)- POS(0))^2+$
(REPLICATE( 1.0, JDIM)#FINDGEN( JDIM)- POS(1))^2)/(2.*SIGMA^2))

return, gauss
closing:
end

;==============================

function mygauss2, idim, jdim, sigmai,sigmaj, pos=pos
 
;-----------------
; parameters check
;-----------------
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, $
   'CALLING SEQUENCE: GAUSS=  mygauss2( idim, jdim, sigmai, sigmaj, pos=pos)'
   GOTO, CLOSING
 ENDIF
if n_elements( pos) eq 0 then begin
  pos= make_array( 2, /float)
  pos[0]= FLOAT(IDIM)/2.0
  pos[1]= FLOAT(JDIM)/2.0
endif
if n_elements( pos) eq 1 then begin
  print, 'Give 2 numbers for POS'
  goto, closing
endif
 
GAUSS= EXP(-(((FINDGEN( IDIM)#REPLICATE( 1.0, IDIM)- POS[0]) ^2   )/(2.*sigmai^2)+$
        (( REPLICATE( 1.0, JDIM)#FINDGEN( JDIM)- POS[1])^2    )/(2.*SIGMAj^2)))
 
return, gauss
closing:
end

;==============================

function mygauss2a, idim, jdim, sigmai,sigmaj, Angle, pos=pos
 
;-----------------
; parameters check
;-----------------
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, $
   'CALLING SEQUENCE: GAUSS=  mygauss2( idim, jdim, sigmai, sigmaj, Angle, pos=pos)'
   GOTO, CLOSING
 ENDIF
if n_elements( pos) eq 0 then begin
  pos= make_array( 2, /float)
  pos[0]= FLOAT(IDIM)/2.0
  pos[1]= FLOAT(JDIM)/2.0
endif
if n_elements( pos) eq 1 then begin
  print, 'Give 2 numbers for POS'
  goto, closing
endif

Y = FINDGEN( IDIM)#REPLICATE( 1.0, IDIM)- POS[0]
X = REPLICATE( 1.0, JDIM)#FINDGEN( JDIM)- POS[1]
A = Angle / 180. * !PI    
X1 =  cos(A )* X  -  sin(A)*Y
Y1 =  sin(A) * X  +  cos(A)*Y
 
GAUSS= EXP(-((X1^2   )/(2.*sigmai^2)+ (Y1^2    )/(2.*SIGMAj^2)))
 
return, gauss
closing:
end

;==============================


pro setpsbw, portrait=portrait, filename=filename
if not keyword_set(filename) then filename = 'idl.ps'

set_plot, 'PS' 
if not keyword_set(portrait) then device, filename=filename, /landscape $
else device, filename=filename, /portrait,/encaps   ; ,xsize=15,ysize=15
end

pro setps, portrait=portrait, filename=filename
if not keyword_set(filename) then filename = 'idl.ps'

set_plot, 'PS' 
if not keyword_set(portrait) then device, filename=filename, /landscape $
else device, filename=filename, /portrait   ; ,xsize=15,ysize=15
end


pro endps
device, /close
set_plot, 'X'
end

;==============================

function origin, frame
  ;; Defines the origin of even- and odd-sized frames.
  n = size(frame)
  result = intarr(2)
  result[0] = n[1]/2
  if (n(0) eq 2) then result[1] = n[2]/2 else result[1] = 0
  return, result
end

;==============================

function centre, image
  ;; Shift the origin of the image from the bottom left corner to the centre.
  o = origin(image)
  s = size(image)
  if (s[0] eq 1) then $
    return, shift(image, o[0]) $
  else $
    return, shift(image, o[0], o[1])
end

;==============================

function decentre, image
  ;; Shift the origin of the image from the centre to the bottom left corner.
  o = origin(image)
  s = size(image)
  if (s[0] eq 1) then $
    return, shift(image, -o[0]) $
  else $
    return, shift(image, -o[0], -o[1])
end

;==============================
function dft, data
  ;; Discrete Fourier transform for centred origin pixel
  ;; The IDL FFT normalizes during the forward transformation.
  ;; return, fft(data, -1)
  return, centre(fft(decentre(data), -1))*N_ELEMENTS(data)
end
;==============================

function dfti, ft
  ;; Inverse discrete Fourier transform for centred origin pixel
  ;; return, fft(ft, +1)
  return, centre(fft(decentre(ft), +1))/N_ELEMENTS(ft)
end
;==============================
function power_spectrum, x
  ;; The squared one
  ft = dft(x)
  return, float(ft*conj(ft))
end

function fftps, x
  ;; The squared one
  ft = dft(x)
  return, float(ft*conj(ft))
end
;==============================
function fftmod, x
  ;; The squared one
  ft = dft(x)
  return, sqrt(float(ft*conj(ft)))
end
;==============================
function fftreal, x
  ;; The squared one
  ft = dft(x)
  return, float(ft)
end
;==============================
function fftimag, x
  ;; The squared one
  ft = dft(x)
  return, float(ft)
end
;==============================
function fftphase, x
  ;; The squared one
  ft = dft(x)
  r = float(ft)
  i = imaginary(ft)
  p = atan(i/r)
  return,p
end

;==============================

function conv, ima, psf
return, float( dfti( dft(ima)*dft(psf) ) )
end

;==============================


function xy, i, array
  ;; Transforms a 1D-index int array into its 2D-equivalent
  n = size(array)
  result = intarr(2)
  result[0] = my_mod(i, n[1])
  result[1] = (i - result[0])/n[1]
  return, result
end

function xyz, i, array
  ;; Transforms a 1D-index int array into its 3D-equivalent
  n = size(array)
  result = intarr(2)
  result[0] = my_mod(i, n[1])
  result[1] = (i - result[0])/n[1]
  return, result
end

;==============================

FUNCTION type_code, object
output=-1
IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=type_code(object)'
   GOTO, CLOSING
 ENDIF
n = size(object)
output = n(n[0]+1)
CLOSING:
RETURN, output
END

FUNCTION type_of, object
output=-1
IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=type_of(object)'
   GOTO, CLOSING
 ENDIF
 tc = type_code(object)
  case tc of
    0: data_type = "UNDEFINED"
    1: data_type = "BYTE"
    2: data_type = "INTEGER"
    3: data_type = "LONGWORD"
    4: data_type = "FLOATING"
    5: data_type = "DOUBLE"
    6: data_type = "COMPLEX"
    7: data_type = "STRING"
    8: data_type = "STRUCTURE"
  endcase
  output = data_type
CLOSING:  
 RETURN, output
 END
 
 pro zzzz
 print, 'zzzz'
 end
 ;==============================
pro make_mr1_help
mk_library_help, "pro", "help/MultiResol.help"
end


pro qpplot, x, y, op=op
print, "OK"
i1 = sort(x)
i2 = sort(y)
vs = size(x)
x1 = x(i1)
y1 = y(i2)
if not keyword_set(op) then plot, x1, y1 else oplot, x1, y1, line=1
; plot, x1, y1, yrange=[-5,5]
end
 ;==============================

function getnsig, eps
   nsig = ABS(inverf(1-eps))*sqrt(2.)
   return, nsig
end

 ;==============================

function geteps, Nsig
  eps = 1. - errorf(double(Nsig) / sqrt(double(2.)))
   return, eps
end

 ;==============================

function strc, x
return, STRCOMPRESS(string(X), /REMOVE_ALL)
end

;==============================

function linspace, A, B, N
L = DINDGEN(N) / (N - 1.0D) * (B - A) + A
return, L
END

;==============================

function logspace, A, B, N
L = DINDGEN(N) / (N - 1.0D) * (B - A) + A
L = 10^L
return, L
END

;==============================

