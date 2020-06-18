pro mk_clean, file=file

; Oct 08 - Written by Anais Rassat
; PURPOSE: Removes all temporary CAMB files.

if not keyword_set(file) then file = './'
command = strcompress('rm '+file+'*cambtemp*')
spawn, command
end
