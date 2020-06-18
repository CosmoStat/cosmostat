
pro mr_nla, ima, tab, ErrTab, tfilter=tfilter, notnla=notnla, write=write, tabx=tabx, nscale=nscale, disp=disp, LC=LC, film=film

if N_PARAMS() LT 1 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr_nla, ima, tab, ErrTab, tfilter=tfilter, notnla=notnla, write=write, tabx=tabx, nscale=nscale, disp=disp, LC=LC, film=film'
  goto, DONE
end

data = ima
data = data /sigma(data)
vs = size(data)
Nx = vs[1]
Ny = vs[2]
ScaleX  = fix(  alog(float(Nx) / 4. * 3.) / alog(2.))  
ScaleY  = fix(  alog(float(Ny) / 4. * 3.) / alog(2.))  

if not keyword_set(nscale) then Nscale = MIN( [ScaleX, ScaleY]) $
else if keyword_set(LC)  then begin
   ScaleX = Nscale
   ScaleY = Nscale
end
lastNx = Nx / 2^(ScaleX-2)
lastNy = Ny / 2^(ScaleY-2)
print, "Nx = ", Nx, " Ny = ", Ny, " NscaleX = ", ScaleX,  " NscaleY = ", ScaleY

if not keyword_set(tabx) then begin
tab=[0.001, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, $
    0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.]
end else begin
tab = tabx
end
n = (size(tab))[1]

;         [-T type_of_filters]
;              1: Biorthogonal 7/9 filters 
;              2: Daubechies filter 4 
;              3: Biorthogonal 2/6 Haar filters 
;              4: Biorthogonal 2/10 Haar filters 
;              5: Odegard 9/7 filters 
;              6: 5/3 filter 
;              7: User's filters 
;              8: Haar filter 
;              9: 3/5 filter 

if  keyword_set(Tfilter) then TFilter=' -T ' + STRCOMPRESS(STRING(Tfilter), /REMOVE_ALL) $
else TFilter=' -T 1'

optn = ' -n ' + STRCOMPRESS(STRING(Nscale), /REMOVE_ALL)
optx = ' -x ' + STRCOMPRESS(STRING(ScaleX), /REMOVE_ALL)
opty = ' -y ' + STRCOMPRESS(STRING(ScaleY), /REMOVE_ALL)

if keyword_set(LC) then OPTLC = TFilter  $
else OPTWT = '-t14 -L ' + TFilter  

ErrTab = fltarr(n)
prem = fltarr(nx,ny,n)
for i=0,n-1 do begin
   if keyword_set(notnla) then begin 
       if keyword_set(LC) then  lc_nfirstinmax, data, rec, tab[i], opt=OPTLC, nsx=ScaleX, nsy=ScaleY  $
       else nfirstinmax, data, rec, tab[i], opt=OPTWT, nscale=Nscale  
   end else begin
        if keyword_set(LC) then nfirstlc, data, rec, nfper=tab[i]*100., opt=OPTLC, nsx=ScaleX, nsy=ScaleY  $
        else nfirstwt, data, rec,  nfper=tab[i]*100., opt=OPTWT, nscale=Nscale
   end
   residu = data-rec
   ErrTab[i] = sigma(residu[10:nx-10,10:ny-10])
   prem[*,*,i] = rec
   if keyword_set(write) then begin
     num = STRCOMPRESS(string(tab[i]), /REMOVE_ALL)
     if keyword_set(LC) then name = 'filmlc_' + num  + '.fits'   $
     else name = 'filmwt_' + num  + '.fits'
     print, name
     writefits,  name, rec
   end
   if keyword_set(disp) then tvscl, rec
   end
   film = prem
   help, film
   if keyword_set(disp) then plot, tab, alog(ErrTab)
   nscale=0
 
DONE:

END

  
