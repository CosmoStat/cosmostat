function df_match_lines, decspectra, estred, rms, input

sz=size(decspectra.data_clean)
ngal = sz[2]

nlinestomatch = n_elements(input.lines)
if strmatch(input.matchtype,'OR') then numcrit = 1
if strmatch(input.matchtype,'AND') then numcrit = nlinestomatch

;; first read in line file

readcol,input.indir+input.linefile, f='d,a', linelam, linename

linepix = round(double(alog10(linelam) - alog10(input.data_lmin))/double(input.lstep))

;; compute the pixel shift corresponding to each galaxy redshift
delpix = round(double(alog10(1.+(estred))/input.lstep))

;; remove galaxies that have no features
excl=where(delpix lt 0, nexcl, complement=incl, ncomplement=nincl)

nmatches = intarr(ngal, nlinestomatch)

;; now loop over galaxies

for i = 0, nincl - 1 do begin
   id = incl[i]

   peaks = df_get_peaks(decspectra.negline[*,id], decspectra.posline[*,id], noise = rms, nsigma = input.nsigmapeak, loc = loc)

   matches = df_linematch(loc, linepix, linename, delpix[id], thresh=input.matchthresh)
   ;help,matches
   ;stop

   nametest = (matches.name)[0]
   if strmatch(strtrim(nametest,2),'0') then nmatches[id,*] = 0 else begin
      for match = 0, n_elements(matches.distance) - 1 do $
         for line = 0, nlinestomatch - 1 do begin
;         print,matches.name[match],input.lines[line]
;         print,strmatch(matches.name[match],input.lines[line])
         if strmatch(matches.name[match],input.lines[line]) then nmatches[id,line] = 1
      endfor
   endelse
   ;print,nmatches[id,*]
   ;read,idum
endfor

nm = total(nmatches,2)

idx = where(nm ge numcrit,count)

clean_catalogue = {indices:idx, redshifts:estred[idx]}

return,clean_catalogue
end

