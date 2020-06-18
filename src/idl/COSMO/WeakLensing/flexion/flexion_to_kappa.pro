pro flexion_to_kappa, inmap, outmap, pixscale, pad = pad, inverse = inverse, $
                      Gflex = Gflex

sz = (size(inmap))
nx = sz[1]
ny = sz[2]
if not keyword_set(inverse) and sz[3] ne 2 then begin
   message, 'Cannot compute kappa from only one flexion component! Returning...',$
            /info
   goto, ret
endif

if not keyword_set(pad) then pad = 1

nxpad = nx*pad
nypad = ny*pad
theta1=nxpad*pixscale ;; pixscale should be in radians
theta2=nypad*pixscale

xm = nxpad/2-nx/2
ym = nypad/2-ny/2

;; pad input data

inmappad = replicate(double(0), nxpad, nypad, 2)
inmappad[((pad-1)*nx/2),((pad-1)*ny/2),*] = inmap

;; compute relevant fourier components

inmappad = fft(inmappad, dimension=[1,2],/overwrite,/double)
l1=indgen(nxpad)#replicate(1,nypad)-(nxpad-1)/2
l1=shift(l1,-(nxpad-1)/2,-(nxpad-1)/2)

l2=replicate(1,nxpad)#indgen(nypad)-(nypad-1)/2
l2=shift(l2,-(nypad-1)/2,-(nypad-1)/2)

l1=2.*!pi/theta1*float(l1)
l2=2.*!pi/theta2*float(l2)

ns1=(pad-1)*nx/2
ns2=(pad-1)*ny/2

inmappad(0,0)=0.


;; compute kappa from flexion

if not keyword_set(inverse) then begin

   if not keyword_set(Gflex) then begin

      kappa = complex(0d,1d)(l1*inmappad[*,*,0]+l2*inmappad[*,*,1])/(l1^2+l2^2)
      
   endif else begin
      
      kappa = complex(0d,1d)*((l1^3-3d*l1*l2^2)*inmappad[*,*,0]+$
                              (l2^3-3*l1^2*l2)*inmappad[*,*,1])/(l1^2+l2^2)^2
      
   endelse
   kappa = (fft(kappa,/double,/inverse,/overwrite))[ns1:ns1+nx-1,$
                                                    ns2:ns2+nx-1]

   outmap = [[[real_part(kappa)]],[[imaginary(kappa)]]]

endif else begin

   if not keyword_set(Gflex) then begin
      
      f1=(double(fft(complex(0d,1d)*l1*inmappad,/double,/inverse)))[ns1:ns1+nx-1,$
                                                                    ns2:ns2+ny-1]
      f2=(double(fft(complex(0d,1d)*l2*inmappad,/double,/inverse)))[ns1:ns1+nx-1,$
                                                                    ns2:ns2+ny-1]
      
      outmap = [[[f1]],[[f2]]]
      
   endif else begin

      g1kern = -2d*complex(0d,1d)*(l1^3-3d*l1*l2^2)/(l1^2+l2^2)
      g2kern = -2d*complex(0d,1d)*(3d*l1^2*l2-l2^3)/(l1^2+l2^2)

      g1 = (double(fft(g1kern*inmappad,/double,/inverse)))[ns1:ns1+nx-1,$
                                                           ns2:ns2+ny-1]
      
      g2 = (double(fft(g2kern*inmappad,/double,/inverse)))[ns1:ns1+nx-1,$
                                                           ns2:ns2+ny-1]

      outmap = [[[g1]],[[g2]]]
   endelse
   

endelse


ret:
return
end
