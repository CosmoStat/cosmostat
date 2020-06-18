PRO pwt1d, signal, trans, nbscale=nbscale, typetrans=typetrans, no_ft=no_ft

signal=reform(signal)
size=size(signal)
npoints=size[1]

;**********************************************************************************

;test puissance de 2

pow2=0
if (ceil(alog(npoints)/alog(2)) eq floor(alog(npoints)/alog(2)) ) then pow2=1
;**********************************************************************************



if not keyword_set(nbscale) then nbscale=4 ; nbscale=4(default) 
if not keyword_set(typetrans) then typetrans=2   ;algo pyramidal decime par defaut


wcoef=dblarr(npoints,nbscale)
imrec=dblarr(npoints,nbscale)


from=fltarr(nbscale)
to=fltarr(nbscale)
wcoef_im=dblarr(npoints,nbscale)


if not keyword_set( no_ft ) then begin	
	
	ligne = fft( signal, -1, /double ) * N_ELEMENTS( signal )
endif else ligne = signal	



CASE 1 OF 

;**********************************************************************************

(typetrans eq 1) and (pow2 eq 1) : BEGIN

	
	;print,''
	;print,''
	
	;print,'              undecimated pyramidal 1D transform'
	;print,'                       nombre de termes en puissance de 2'
	
	for numscale=0, nbscale-2 do begin

		
		w = dcomplexarr( npoints, nbscale )
		t = dcomplexarr( npoints, nbscale )
		
		compute_filter_h, npoints, 1, h   ;passe bas
		compute_filter_g, npoints, 1, g   ;passe haut
		
		w= ligne * g
		wcoef( *, numscale ) = real_part( fft( w, 1, /double ) / N_ELEMENTS(w) )     ;gives the wavelet coefficients at scale 2^j 
		
		t= ligne * h
		imrec( *, numscale ) = real_part( fft( t, 1, /double ) / N_ELEMENTS(t) ) ;gives the image at scale 2^(j+1)
		 
		from( numscale ) = 0
		to( numscale ) = npoints
		
		ligne = dcomplexarr( npoints )
		ligne = fft( imrec( *, numscale), -1, /double) * N_ELEMENTS( imrec( *,numscale ))
		
	
	endfor
	
	
	to( nbscale-1 ) = to( nbscale-2 )
	

END



;**********************************************************************************


(typetrans eq 1) and (pow2 eq 0) : BEGIN

	
	;print,''
	;print,'     Transformee Ondelettes 1D Non Decimee'
	

	for numscale = 0, nbscale-2 do begin

		
		w = dcomplexarr( npoints, nbscale )
		t = dcomplexarr( npoints, nbscale )
		
		
		compute_filter_h, npoints, 1, h   ;passe bas
		compute_filter_g, npoints, 1, g   ;passe haut
		
		w = ligne * g
		wcoef(*,numscale) = real_part( fft( w,1, /double) / N_ELEMENTS(w))     ;gives the wavelet coefficients at scale 2^j 
		
		t = ligne * h
		
	
		imrec(*,numscale) = real_part( fft( t, 1, /double) / N_ELEMENTS(t)) ;gives the image at scale 2^(j+1)
		
		from(numscale) = 0
		to(numscale) = npoints
		
		ligne = dcomplexarr( npoints / 2 )
		ligne = fft( imrec(*, numscale), -1, /double ) * N_ELEMENTS( imrec(*,numscale))
		
	
	endfor
	
	to(nbscale-1) = to(nbscale-2)

	
END



;**********************************************************************************



(typetrans eq 2) and (pow2 eq 1) : BEGIN


	;print,''
	;print,''
	;print,'                  Transformee Ondelettes methode pyramidale decimee '
	;print,'                  (Le nombre de termes est en puissance de 2)'
	
	for numscale = 0, nbscale-2 do begin

		
		w = dcomplexarr( npoints, nbscale )
		t = dcomplexarr( npoints, nbscale )
		tdecim = dcomplexarr( npoints/2 )
		
		compute_filter_h, npoints, 1, h   ;passe bas
		compute_filter_g, npoints, 1, g   ;passe haut
		
		w = ligne * g
		wcoef( 0:npoints-1, numscale) = real_part( fft( w, 1, /double ) / N_ELEMENTS(w))     ;gives the wavelet coefficients at scale 2^j 
		


		t = ligne * h
		ind=[ findgen( npoints/4 ), findgen( npoints/4 ) + 3 * npoints/4]
		tdecim=t(ind)
		imrec( 0:npoints/2-1, numscale ) = real_part( fft( tdecim, 1, /double ) / ( 2 * N_ELEMENTS(tdecim))) ;gives the image at scale 2^(j+1)
		
		from(numscale) = 0
		to(numscale) = npoints
		
		
		ligne = dcomplexarr( npoints / 2 )
		ligne = fft( imrec( 0:npoints/2-1, numscale ), -1, /double) * npoints / 2
		npoints = npoints / 2
	
	endfor
	
	to( nbscale - 1 ) = to( nbscale - 2 )/ 2

END

;**********************************************************************************


(typetrans eq 2) and (pow2 eq 0) : BEGIN


	;print,''
	;print,''
	;print,'          Transformee Ondelettes methode pyramidale decimee '
	
	
	for numscale = 0, nbscale-2 do begin

		
		w = dcomplexarr( npoints, nbscale )
		t = dcomplexarr( npoints, nbscale )
		
		tdecim = dcomplexarr( ceil( npoints / 2.) )
		
		compute_filter_h, npoints, 1, h   ;passe bas
		compute_filter_g, npoints, 1, g   ;passe haut
		
		w = ligne * g
		wcoef( 0:npoints-1, numscale ) = real_part( fft( w, 1, /double ) / N_ELEMENTS(w) )     ;gives the wavelet coefficients at scale 2^j 
		
		t = ligne * h
		ind = [findgen(floor(ceil(npoints/2.)/2)),findgen(ceil(npoints/4.))+npoints-ceil(npoints/4.)]  
		tdecim = t(ind)
		imrec( 0:ceil(npoints/2.)-1, numscale ) = real_part( fft(tdecim,1,/double) / ( 2*N_ELEMENTS(tdecim)))
		
		from(numscale) = 0
		to(numscale) = npoints
		
		ligne = dcomplexarr( ceil( npoints / 2. ))
		ligne = fft(imrec(0:ceil( npoints / 2.)-1,numscale),-1,/double)*(ceil(npoints/2.))
		npoints = ceil(npoints/2.)
	
	
	endfor
	
	to(nbscale-1) = ceil(to(nbscale-2)/2.)

END

ENDCASE
;**********************************************************************************


for numscale = 0, nbscale - 2 do begin
	
	wcoef_im( *, numscale ) = wcoef( *, numscale )
endfor



wcoef_im( *,nbscale-1 ) = imrec( *, nbscale - 2 )
trans={ n_band : nbscale, coef : wcoef_im, typetrans : typetrans, from : from,to : to}



;print,''
END
