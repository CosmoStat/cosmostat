
pro test_anomalies
;.com Stat_Anomalies
; Mai 2013
; Test la version C++/IDL du code StatAnomalies

;Example of how to use code: 
t = getcmb() & t = t/1d3 ; so that entry is in mK
fidstat = set_fid_stat()
Results = Stat_Anomalies(t, fidstat=fid)
stop

; Below you can see a more detailed example: 

nside = 512 ; nside you want to resize maps to
nsiderot=2 ; nside of rotations for tests

t0 =rims('../WMAP/toh_wmap_cleaned_yr1_v1_NESTEDv2.fits')
t0 = reform(t0[*,0])
t0 = reform(mrs_resize(t0,nside=nside))

t= rims('../WMAP/wmap_ilc_3yr_v2.fits')
t = reform(t[*,0])
t = reform(mrs_resize(t,nside=nside))


t2= rims('../WMAP/wmap_ilc_5yr_v3.fits')
t2 = reform(t2[*,0])
t2 = reform(mrs_resize(t2,nside=nside))

t3 = rims('../WMAP/wmap_ilc_7yr_v4.fits')
t3 = reform(t3[*,0])
t3 = reform(mrs_resize(t3,nside=nside))

t4 = rims('../WMAP/wmap_ilc_9yr_v5.fits')
t4 = reform(t4[*,0])
t4 = reform(mrs_resize(t4,nside=nside))

name = ['W1','W3', 'W5', 'W7', 'W9']
map = dblarr(n_elements(t), n_elements(name))
map[*,0] = t0 & map[*,1] =t & map[*,2] =t2 & map[*,3] =t3 & map[*,4]=t4 ; maps should be in units of microK
data = {name:name, map:map}
help, data,/st

; Set up the fiducial structure with all options for calculations of
; statistical anomalies
; If you don't know what to input you can just do: 
; >fidstat = set_fid_stat() ; with no explicit options
fidstat = set_fid_stat(nsiderot=nsiderot, lowquad={test:1}, quadoct={test:1}, planaroct={test:1}, AOE={test:1}, Mparity={test:1} ) 

print, 'Starting to calculate rotations'
for i = 0, n_elements(data.name)-1 do begin
   a = Stat_Anomalies(reform(data.map[*,i]), fidstat=fidstat,mapname=data.name[i]) 
   if i eq 0 then Results = replicate( a, n_elements(data.name))
   Results[i] = a
endfor

stop

end
