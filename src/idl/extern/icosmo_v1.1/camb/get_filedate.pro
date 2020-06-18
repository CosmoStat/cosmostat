pro get_filedate,filename

; June 08 - Written by Anais Rassat
; PURPOSE: Returns a string with the year/month/day/hour/minute/second

s=bin_date()
filename = strcompress(string(s(0))+string(s(1))+string(s(2))+string(s(3))+string(s(4))+string(s(5))+'_cambtemp',/remove_all)
print, filename
end
