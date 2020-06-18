pro psnr, ima1, ima2

resi = ima1 - ima2
n = N_ELEMENTS(resi)

print, "sigma(Err) = ", sigma(Resi)


V1 = total(resi^2)/ N

print, "PSNR_dB = ", 10. * alog10( 255.^2 / V1)
print, "SNR_dB = ", 10. * alog10( total(ima1^2) / total(resi^2))

end
