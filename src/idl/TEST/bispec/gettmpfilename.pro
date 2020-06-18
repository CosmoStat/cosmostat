
function gettmpfilename
; spawn,'mktemp',filename_tmp
if  !version.os EQ 'linux' then spawn,'mktemp -p /dev/shm', filename_tmp $
else  spawn,'mktemp',filename_tmp 

; filename_tmp = strcompress('tmp_in_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
return, filename_tmp
end
