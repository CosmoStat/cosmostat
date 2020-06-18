
pro pg, genus, title=title, dim=dim

if keyword_set(dim) then N = float(DIM)^3 else N = 1
N = 256.^3
g = -genus(*,3)
g = g / max(g)
xrang=[-2.,2.]
M=1.5
yrang=[-0.75,1.25]
plot, genus(*,6), g, yticks=5, xthick=3,xrange=xrang, yrange=yrang,  $
         ticklen=0.05,yminor=2,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
         thick=1, xticklen=0.1, ythick=2,/ynozero, title=title
end


pro pg4, g1,g2,g3,g4, title=title, tablegend=tablegend,dx=dx 
g = -g1(*,3)
g = g / max(abs(g))
xrang=[-2.,2.]
yrang=[-0.5,1.]
plot, g1(*,6), g, yticks=4, xrange=xrang, yrange=yrang,  $
         ticklen=0.04,yminor=3,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
         thick=1, xticklen=0.04, xthick=2, ythick=2,/ynozero, title=title, $
	 xtitle='Nu', ytitle='Normalized Genus'

g = -g2(*,3)
g = g / max(g)
oplot, g2(*,6), g, line=1

g = -g3(*,3)
g = g / max(g)
oplot, g3(*,6), g, line=2

g = -g4(*,3)
g = g / max(g)
oplot, g4(*,6), g, line=3

posy=1.08
posx=0.6
delta = 0.13
deltax = 0.27
if keyword_set(dx) then posx=posx+dx
Tab=['Scale 1', 'Scale 2', 'Scale 3', 'Scale 4']
if not keyword_set(tablegend) then tablegend=tab
for b=0,3 do   begin
      y = b*delta+posy
      oplot, [posx,posx+deltax], [posy,posy], line=b
      xyouts, posx+deltax+0.1, posy, tablegend[b], charsize=1
      posy = posy - delta
    end
end

pro pg3, g1,g2,g3,title=title, tablegend=tablegend,dx=dx 
g = -g1(*,3)
g = g / max(abs(g))
xrang=[-2.,2.]
yrang=[-0.5,1.]
plot, g1(*,6), g, yticks=4, xrange=xrang, yrange=yrang,  $
         ticklen=0.04,yminor=3,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
         thick=1, xticklen=0.04, xthick=2, ythick=2,/ynozero, title=title, $
	 xtitle='Nu', ytitle='Normalized Genus'

g = -g2(*,3)
g = g / max(g)
oplot, g2(*,6), g, line=1

g = -g3(*,3)
g = g / max(g)
oplot, g3(*,6), g, line=2


posy=1.08
posx=0.6
delta = 0.13
deltax = 0.27
if keyword_set(dx) then posx=posx+dx
Tab=['Scale 1', 'Scale 2', 'Scale 3']
if not keyword_set(tablegend) then tablegend=tab
for b=0,2 do   begin
      y = b*delta+posy
      oplot, [posx,posx+deltax], [posy,posy], line=b
      xyouts, posx+deltax+0.1, posy, tablegend[b], charsize=1
      posy = posy - delta
    end
end


;==========================================================================
pro makefig, ps=ps
n0 = rim('g_n1_c0.fits')                   
n1 = rim('g_n1_c1.fits')                
n2 = rim('g_n1_c2.fits')                  
n3 = rim('g_n1_c4.fits')                  
n4 = rim('g_n1_c8.fits')

if keyword_set(PS) then setps, filename='fig_genus_grfn1.ps', /portrait
pg4, n1, n2,n3, n4, title='Gaussian Random Field (n=-1)', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wn1 = rim('g_wt_n1_band_1.fits')                
wn2 = rim('g_wt_n1_band_2.fits')                  
wn3 = rim('g_wt_n1_band_3.fits')                  
wn4 = rim('g_wt_n1_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_grfn1.ps', /portrait
pg4, wn1, wn2, wn3, wn4, title='Gaussian Random Field (n=-1)', $
     tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps

v1 = rim('g_virgo128_c1.fits')
v2 = rim('g_virgo128_c2.fits')
v3 = rim('g_virgo128_c4.fits')
v4 = rim('g_virgo128_c8.fits')
if keyword_set(PS) then setps, filename='fig_genus_virgo.ps', /portrait
pg4, v1, v2,v3, v4, title='Virgo Simulation', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wv1 = rim('g_wt_virgo128_band_1.fits')
wv2 = rim('g_wt_virgo128_band_2.fits')
wv3 = rim('g_wt_virgo128_band_3.fits')
wv4 = rim('g_wt_virgo128_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_virgo.ps', /portrait
pg4, wv1, wv2, wv3, wv4, title='Virgo Simulation', $
     tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps

g1 = rim('g_gif_c1.fits')                
g2 = rim('g_gif_c2.fits')                  
g3 = rim('g_gif_c4.fits')                  
g4 = rim('g_gif_c8.fits')
if keyword_set(PS) then setps, filename='fig_genus_gif.ps', /portrait
pg4, g1, g2,g3, g4, title='GIF Simulation', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wg1 = rim('g_wt_gif_band_1.fits')                
wg2 = rim('g_wt_gif_band_2.fits')                  
wg3 = rim('g_wt_gif_band_3.fits')                  
wg4 = rim('g_wt_gif_band_4.fits')
if keyword_set(PS) then setps, filename='fig_genus_wave_gif.ps', /portrait
pg4, wg1, wg2,wg3, wg4, title='GIF Simulation', $
      tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps

p1 = rim('g_poisson_c1.fits')                
p2 = rim('g_poisson_c2.fits')                  
p3 = rim('g_poisson_c4.fits')                  
p4 = rim('g_poisson_c8.fits')
if keyword_set(PS) then setps, filename='fig_genus_poisson.ps', /portrait
pg4, p1, p2,p3, p4, title='Poisson Simulation', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wp1 = rim('g_wt_poisson_band_1.fits')                
wp2 = rim('g_wt_poisson_band_2.fits')                  
wp3 = rim('g_wt_poisson_band_3.fits')                  
wp4 = rim('g_wt_poisson_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_poisson.ps', /portrait
pg4, wp1, wp2,wp3, wp4, title='Poisson Simulation', $
      tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps
end
     

;===============================================================

pro makefig2, ps=ps

n1 = rim('g_clus_c1.fits')                
n2 = rim('g_clus_c2.fits')                  
n3 = rim('g_clus_c4.fits')                  
n4 = rim('g_clus_c8.fits')

if keyword_set(PS) then setps, filename='fig_genus_clus.ps', /portrait
pg4, n1, n2,n3, n4, title='clusters', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wn1 = rim('g_wt_clus_band_1.fits')                
wn2 = rim('g_wt_clus_band_2.fits')                  
wn3 = rim('g_wt_clus_band_3.fits')                  
wn4 = rim('g_wt_clus_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_clus.ps', /portrait
pg4, wn1, wn2, wn3, wn4, title='clusters', $
     tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps


n1 = rim('g_fil_c1.fits')                
n2 = rim('g_fil_c2.fits')                  
n3 = rim('g_fil_c4.fits')                  
n4 = rim('g_fil_c8.fits')

if keyword_set(PS) then setps, filename='fig_genus_fil.ps', /portrait
pg4, n1, n2,n3, n4, title='filaments', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wn1 = rim('g_wt_fil_band_1.fits')                
wn2 = rim('g_wt_fil_band_2.fits')                  
wn3 = rim('g_wt_fil_band_3.fits')                  
wn4 = rim('g_wt_fil_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_fil.ps', /portrait
pg4, wn1, wn2, wn3, wn4, title='filaments', $
     tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps

n1 = rim('g_wall_c1.fits')                
n2 = rim('g_wall_c2.fits')                  
n3 = rim('g_wall_c4.fits')                  
n4 = rim('g_wall_c8.fits')

if keyword_set(PS) then setps, filename='fig_genus_wall.ps', /portrait
pg4, n1, n2,n3, n4, title='walls', $
     tablegend=['Sigma=1','Sigma=2', 'Sigma=4','Sigma=8']
if keyword_set(PS) then endps

wn1 = rim('g_wt_wall_band_1.fits')                
wn2 = rim('g_wt_wall_band_2.fits')                  
wn3 = rim('g_wt_wall_band_3.fits')                  
wn4 = rim('g_wt_wall_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_wave_wall.ps', /portrait
pg4, wn1, wn2, wn3, wn4, title='walls', $
     tablegend=['Scale 1','Scale 2', 'Scale 3','Scale 4']
if keyword_set(PS) then endps
end


;===============================================================

pro makefig3, ps=ps

c1 = rim('g_clus_c1.fits')                
c2 = rim('g_clus_c2.fits')                  
c3 = rim('g_clus_c4.fits')                  
c4 = rim('g_clus_c8.fits')
f1 = rim('g_fil_c1.fits')                
f2 = rim('g_fil_c2.fits')                  
f3 = rim('g_fil_c4.fits')                  
f4 = rim('g_fil_c8.fits')
w1 = rim('g_wall_c1.fits')                
w2 = rim('g_wall_c2.fits')                  
w3 = rim('g_wall_c4.fits')                  
w4 = rim('g_wall_c8.fits')

wtc1 = rim('g_wt_clus_band_1.fits')                
wtc2 = rim('g_wt_clus_band_2.fits')                  
wtc3 = rim('g_wt_clus_band_3.fits')                  
wtc4 = rim('g_wt_clus_band_4.fits')
wtf1 = rim('g_wt_fil_band_1.fits')                
wtf2 = rim('g_wt_fil_band_2.fits')                  
wtf3 = rim('g_wt_fil_band_3.fits')                  
wtf4 = rim('g_wt_fil_band_4.fits')
wtw1 = rim('g_wt_wall_band_1.fits')                
wtw2 = rim('g_wt_wall_band_2.fits')                  
wtw3 = rim('g_wt_wall_band_3.fits')                  
wtw4 = rim('g_wt_wall_band_4.fits')

if keyword_set(PS) then setps, filename='fig_genus_3_sig1.ps', /portrait
pg3, c1, f1, w1, title='Genus sigma=1', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_sig2.ps', /portrait
pg3, c2, f2, w2, title='Genus sigma=2', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_sig4.ps', /portrait
pg3, c3, f3, w3, title='Genus sigma=4', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_sig8.ps', /portrait
pg3, c4, f4, w4, title='Genus sigma=8', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps


if keyword_set(PS) then setps, filename='fig_genus_3_wt1.ps', /portrait
pg3, wtc1, wtf1, wtw1, title='Genus Scale=1', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_wt2.ps', /portrait
pg3, wtc2, wtf2, wtw2, title='Genus Scale=2', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_wt3.ps', /portrait
pg3, wtc3, wtf3, wtw3, title='Genus Scale=3', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps
if keyword_set(PS) then setps, filename='fig_genus_3_wt4.ps', /portrait
pg3, wtc4, wtf4, wtw4, title='Genus Scale=4', $
     tablegend=['cluster','filament', 'wall']
if keyword_set(PS) then endps

end
