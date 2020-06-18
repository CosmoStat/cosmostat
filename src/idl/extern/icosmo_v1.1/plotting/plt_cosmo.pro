pro plt_cosmo, cosmo, xvariable, yvariable, overplot=overplot, color=color,_extra=extra

; Oct 2008 - Modified by An.R. Corrected ytitle output.
; Aug. 2008 - Modified by Anais Rassat. Plot volume.
; July 2008 - Modified by Anais Rassat. Color keyword added and
;            _strict_extra changed to _extra as the former does not
;            always work
; July 2008 - Written by Anais Rassat
; ***************// HELP BEGIN //**************
; PURPOSE: Plots any cosmological quantity in cosmo.evol, outputted by
;          mk_cosmo (version 0.10.0 and above)
; INPUTS:  cosmo -   strucuture output by mk_cosmo (new version)
;          xvariable - string indicating xvariable. Accepts 'z' or 'a'
;          yvariable - string indicating yvaraible. Accepts:
;                      'hc', 'chi', 'sk', 'da', 'dl', 
;                      'd','dlnddlna', 'w_a','omega_a',
;                      'omega_m_a', 'omega_l_a', 'omega_k_a', 'dzdr', 'vol'
; OPTIONAL INPUTS:
;          overplot - if selected, does not create a new plot but
;                     overplots on previous. Means the xyvariable
;                     should be the SAME as in the previous plot.
;          color - keyword to change color of line (not axis)
;          _extra - any keyword accepted by plot. If extra keywords are
;                   not accepted by plot, then these are 'quietly ignored'. 
; OUTPUT: None
; OPTIONAL OUTPUT: None
;
; Example 1:
;       > plt_cosmo, cosmo, 'z', 'omega_m_a'
; Example 2:
;       > plt_cosmo, cosmo, 'z', 'omega_l_a', /over, linestyle=3
; ***************// HELP END //**************
tek_color

;Find out which variable wants to be plotted
void = tag_exist(cosmo.evol, xvariable, index=indx)
void = tag_exist(cosmo.evol, yvariable, index=indy)
;Defining what to plot
x = cosmo.evol.(indx)
y = cosmo.evol.(indy)
;Correcting for units for chi and sk:
if yvariable eq 'chi' or yvariable eq 'sk' then y = cosmo.evol.(indy)*cosmo.const.r0

;Define the title on the y axis
CASE yvariable OF
   'hc':  ytitle='H(z) [km.s!U-1!N.Mpc!U-1!N]'
   'chi': ytitle='!7v!X(z) [!NMpc]'
   'sk':  ytitle ='r(z) [!NMpc]'
   'da':  ytitle ='D!DA!N(z) [Mpc]'
   'dl':  ytitle='D!Dl!N(z) [Mpc]'
   'd' :  ytitle='Growth D(z)'
   'dlnddlna': ytitle='Growth Rate dlnD/dlna (z)'
   'w_a': ytitle='w(z)'
   'omega_a': ytitle='!7X!X!Dm!N(z) + !7X!X!D!7K!X!N(z)'
   'omega_m_a':ytitle='!7X!X!Dm!N(z)'
   'omega_l_a':ytitle='!7X!X!D!7K!X!N(z)'
   'omega_k_a':ytitle='!7X!X!Dk!N(z)'
   'dzdr': ytitle='dz/dr (z)'
   'vol': ytitle= ' Volume(z) [Mpc^3]'
;In case people want a or z as the yvariable
   'a': ytitle='a(z)'
   'z': ytitle='z(a)'
ENDCASE

if keyword_set(overplot) then oplot, x, y, color=color,_extra=extra else begin 
   plot, x, y, xtitle=xvariable, ytitle=ytitle, _extra=extra,/nodata
   oplot, x, y, color=color, _extra=extra
endelse

   
   
End
