function read_wigglez_data, filename

;; filename should specify full filepath to data

inspec = mrdfits(filename, 0, header0)
innoise = mrdfits(filename, 1, header1)

refpix = sxpar(header0,'CRPIX1')
refval = sxpar(header0,'CRVAL1')
delta = sxpar(header0,'CDELT1')

;print,refpix, refval, delta

lambda = (findgen(n_elements(inspec))-refpix)*delta+refval
;print,min(lambda), max(lambda)
;read,idum

output = {spectrum:inspec, variance:innoise, lambda:lambda, refpix:refpix, $
          refval:refval, delta:delta, header0:header0, header1:header1}

return, output
end
