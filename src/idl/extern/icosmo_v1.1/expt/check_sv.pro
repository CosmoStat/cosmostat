pro check_sv,sv, probe

;Sep 2008 - Written by Anais Rassat
;PURPOSE: checks if a given survey was designed for a given probe

if probe eq 'sne' then p = 'Supernovae'
if probe eq 'bao' then p = 'BAO'
if probe eq 'lens' then p = 'Weak Lensing'

if where(strmatch(sv.probes, probe)  eq 1) eq -1 then begin
      print, strcompress('check_sv: Survey not defined for '+p+' calculation.')
      print, 'check_sv: Programme stopped.'
stop
endif

end
