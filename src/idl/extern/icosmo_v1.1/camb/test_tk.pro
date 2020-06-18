pro test_tk

f = set_fiducial()
cosmo = mk_cosmo(f)

wt_params_camb_general, f, fileparams=file, filetk=filetk, filepk=filepk
;command =strcompress('/Users/arassat/Work/CAMB/./camb '+string(f))
f.calc.cambpath ='~/Work/CAMB_0902/CAMB'
command =strcompress(f.calc.cambpath+'/./camb', /remove_all)
command = strcompress(command+' '+string(file))
spawn, command


readcol,filetk, k, a, b, c, d, e, tk
readcol, filepk, k2, pk

tk2 = interpol(tk,k,k2,/spline)

plot, k2,pk,/xlo,/ylo,xrange=[1d-3,1d2], yrange=[1d-4,1d5], charsize=2 
oplot, k2, tk2/max(tk2), linestyle=1
oplot, k2, pk/tk2^2*max(tk2)^2/k2
;oplot, 

stop
end
