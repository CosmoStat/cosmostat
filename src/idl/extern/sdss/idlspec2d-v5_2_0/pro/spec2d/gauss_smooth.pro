
function gauss_smooth, x, y, sigma, xout
;
;	Input data in x and y
;	smooth with sigma (in units of x) and evaluate at xout
;
	var = sigma*sigma
	denom = 1.0/sqrt(2.0*!Pi * var)

	yout = xout * 0.0
	
	for i=0,n_elements(xout)-1 do begin
	  ;
	  ;	Find all x's within 10 sigma
	  ;
	  diff = (x*1.d0-xout[i])^2
	  thesex = where(diff LT 100.0*var)
	  if (thesex[0] NE -1) then begin
	    kernel = exp(-0.5*diff[thesex]/var)*denom
	    yout[i] = total(kernel*y[thesex])/total(kernel)
	  endif
	endfor

	return, yout
end

