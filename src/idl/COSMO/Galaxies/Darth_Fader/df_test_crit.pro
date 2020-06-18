pro df_test_crit, fail_only=fail_only, success_only=success_only, thresh=thresh

indata = readfits('./df_data/components.fits')
rms = readfits('./df_data/errorcurve.fits.gz')
lstep = double(alog10(6003d) - alog10(6000d))
lmin = 3600.4782 
sz=size(indata)
nlam = sz[1]

logwl = double(alog10(lmin) + dindgen(nlam)*lstep)

denoised_lines = indata[*,*,0] + indata[*,*,1]

decspectra = {posline:indata[*,*,0], negline:indata[*,*,1]}

redest = readfits('./df_data/redshifts_gen2.fits')

truered = readfits('./df_data/trueredshifts.fits.gz')
readcol,'lines.txt', f='d,a', linelam, linename

linepix = round(double(alog10(linelam) - alog10(lmin))/double(lstep))

;; compute the pixel shift corresponding to each galaxy redshift
delpix = round(double(alog10(1.+(redest))/lstep))
deltrue = round(double(alog10(1.+(truered))/lstep))

excl=where(delpix eq max(delpix), nexcl, complement=incl, ncomplement=nincl)

x=dblarr(nlam)*0
ftab = x
ftab[round(double(alog10(4000d) - alog10(lmin))/double(lstep))] = 1.

fail = intarr(sz[2])
;success = fail
fail[where(abs(redest-truered) gt 0.005*(1+truered))] = 1
success = 1-fail

nmatches = intarr(n_elements(truered))
hdline = nmatches
hgline = nmatches
hbline = nmatches
haline = nmatches
o2line = nmatches
o31line = nmatches
o32line = nmatches
breakline = nmatches
s21line = nmatches
s22line = nmatches

for i = 0, nincl - 1 do begin

;   if keyword_set(fail_only) then if success[incl[i]] then goto, skip
;   if keyword_set(success_only) then if fail[incl[i]] or $
;   (peaks.allpeaks)[incl[i]] ge 6 then goto, skip

   peaks = df_get_peaks(indata[*,incl[i],1], indata[*,incl[i],0], noise=rms, nsigma=0.01, loc=loc)

   matches = df_linematch(loc, linepix, linename, delpix[incl[i]],thresh=thresh)
;   stop

;   print,i
   x=x*0.
   y=x
   dum = where(linepix + delpix[incl[i]] ge 0)
   x[linepix[dum]+delpix[incl[i]]] = 2.*max(abs(denoised_lines[*,incl[i]]))
   y[linepix[dum]+deltrue[incl[i]]] = -2.*max(abs(denoised_lines[*,incl[i]]))

;   loadct,5
;   plot,denoised_lines[*,incl[i]], yrange=[-max(abs(denoised_lines[*,incl[i]])), max(abs(denoised_lines[*,incl[i]]))]
;   oplot,x,color=100
;   oplot,y,color=100,linestyle=1
;   print,"z_est = ", redest[incl[i]], " z_true = ", truered[incl[i]]
;   print,"Peak count = ", (peaks.allpeaks)
;   if success[incl[i]] then print, 'Correct' else print, 'Failed'
;   oplot,shift(ftab,delpix[incl[i]])*2*max(abs(denoised_lines[*,incl[i]])),color=200
;   oplot,-shift(ftab,deltrue[incl[i]])*2*max(abs(denoised_lines[*,incl[i]])),color=200,linestyle=1
   
;   print,"lines matched, distance in pixels"
;   print,matches.name, matches.distance

   nametest=(matches.name)[0]
   if strmatch(string(nametest),'0') then nmatches[incl[i]] = 0 else begin
      nmatches[incl[i]] = n_elements(matches.distance)
      for match = 0, n_elements(matches.distance)-1 do begin
         if strmatch((matches.name)[match],"H_delta") then hdline[incl[i]] = 1
         if strmatch((matches.name)[match], "H_gamma") then hgline[incl[i]] = 1 
         if strmatch((matches.name)[match],"H_beta") then hbline[incl[i]] = 1
         if strmatch((matches.name)[match],"Break") then breakline[incl[i]] = 1
         if strmatch((matches.name)[match],"OII_1") then o2line[incl[i]] = 1
         if strmatch((matches.name)[match],"OIII_1") then o31line[incl[i]] = 1
         if strmatch((matches.name)[match],"OIII_2") then o32line[incl[i]] = 1
         if strmatch((matches.name)[match],"SII_1") then s21line[incl[i]] = 1
         if strmatch((matches.name)[match],"SII_2") then s22line[incl[i]] = 1
      endfor
   endelse
;   if idum eq 2 then goto,skip2
skip: 
endfor

print,'Analysing results'
print,'N_successes = ', n_elements(where(success))

idx=where(o2line gt 0 or breakline gt 0, count)

print,'Retaining on OII or Break'
print,'Retention rate', double(count)/double(2860)*100d
print,'Success rate after cleaning', double(n_elements(where(success[idx])))/double(count)*100d
stop

print,'Results in single peak cases'
idx=where(nmatches eq 1,count)
print,'Number, successes, failures', count, n_elements(where(success[idx])), n_elements(where(fail[idx]))
hline = where(haline[idx] gt 0 or hbline[idx] gt 0 or hgline gt 0 or hdline gt 0,count)
print,'Hydrogen line, successes, failures', count, n_elements(where(success[idx[hline]])), n_elements(where(fail[idx[hline]]))
bline = where(breakline[idx] gt 0,count)
print, 'Break, successes, failures', count,n_elements(where(success[idx[bline]])), n_elements(where(fail[idx[bline]]))
o2l = where(o2line[idx] gt 0, count)
print,'OII line, successes, failures', count, n_elements(where(success[idx[o2l]])), n_elements(where(fail[idx[o2l]]))
o3l = where(o31line[idx] gt 0 or o32line[idx],count)
print,'OIII line, successes, failures', count, n_elements(where(success[idx[o3l]])), n_elements(where(fail[idx[o3l]]))

print,'Results in 2 peak cases'
idx=where(nmatches eq 2,count)
print,'Number, successes, failures', count, n_elements(where(success[idx])), n_elements(where(fail[idx]))
hline = where(haline[idx] gt 0 or hbline[idx] gt 0 or hgline gt 0 or hdline gt 0,count)
print,'Hydrogen line, successes, failures', count, n_elements(where(success[idx[hline]])), n_elements(where(fail[idx[hline]]))
bline = where(breakline[idx] gt 0,count)
print, 'Break, successes, failures', count,n_elements(where(success[idx[bline]])), n_elements(where(fail[idx[bline]]))
o2l = where(o2line[idx] gt 0, count)
print,'OII line, successes, failures', count, n_elements(where(success[idx[o2l]])), n_elements(where(fail[idx[o2l]]))
o3l = where(o31line[idx] gt 0 or o32line[idx],count)
print,'OIII line, successes, failures', count, n_elements(where(success[idx[o3l]])), n_elements(where(fail[idx[o3l]]))

print,'Results in 3 peak cases'
idx=where(nmatches eq 3,count)
print,'Number, successes, failures', count, n_elements(where(success[idx])), n_elements(where(fail[idx]))
hline = where(haline[idx] gt 0 or hbline[idx] gt 0 or hgline gt 0 or hdline gt 0,count)
print,'Hydrogen line, successes, failures', count, n_elements(where(success[idx[hline]])), n_elements(where(fail[idx[hline]]))
bline = where(breakline[idx] gt 0,count)
print, 'Break, successes, failures', count,n_elements(where(success[idx[bline]])), n_elements(where(fail[idx[bline]]))
o2l = where(o2line[idx] gt 0, count)
print,'OII line, successes, failures', count, n_elements(where(success[idx[o2l]])), n_elements(where(fail[idx[o2l]]))
o3l = where(o31line[idx] gt 0 or o32line[idx],count)
print,'OIII line, successes, failures', count, n_elements(where(success[idx[o3l]])), n_elements(where(fail[idx[o3l]]))

print,'Results in 4 peak cases'
idx=where(nmatches eq 4,count)
print,'Number, successes, failures', count, n_elements(where(success[idx])), n_elements(where(fail[idx]))
hline = where(haline[idx] gt 0 or hbline[idx] gt 0 or hgline gt 0 or hdline gt 0,count)
print,'Hydrogen line, successes, failures', count, n_elements(where(success[idx[hline]])), n_elements(where(fail[idx[hline]]))
bline = where(breakline[idx] gt 0,count)
print, 'Break, successes, failures', count,n_elements(where(success[idx[bline]])), n_elements(where(fail[idx[bline]]))
o2l = where(o2line[idx] gt 0, count)
print,'OII line, successes, failures', count, n_elements(where(success[idx[o2l]])), n_elements(where(fail[idx[o2l]]))
o3l = where(o31line[idx] gt 0 or o32line[idx],count)
print,'OIII line, successes, failures', count, n_elements(where(success[idx[o3l]])), n_elements(where(fail[idx[o3l]]))


stop
skip2:
end
