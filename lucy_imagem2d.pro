pro lucy_imagem2d
pasta = 'C:\Users\incta\Desktop\share\5643_cubes\calib_5643'
arquivo = 'bulge.fits'
image = mrdfits(pasta+'\'+arquivo)
; Coloque a FWHM em arcsec. Se tiver apenas em px, coloque a pixelagem em 1
;pixelagem = 1.0; original 0.076 nicmos e 0.045 WFPC
;fwhm = 2.9
niter = 6
z = size(image)

print, z
;psf = psf_gaussian(NPIXEL=[z[1],z[2]], FWHM = fwhm / pixelagem, /normalize)
;psf = psf_airy_v3(z[1], z[2], 6001.281E-10, pixelagem, 2.4)
psf = mrdfits('C:\Users\incta\Desktop\share\5643_cubes\calib_5643\psf_final_4std.fits')
for k = 1, niter do Max_Likelihood, image[*,*,0], psf, deconvlike
itstrg = string(niter)
itstrg=STRTRIM(itstrg,1)
;stop


mwrfits, deconvlike, pasta+'\'+'deconvolved1'+itstrg+arquivo+'.fits'
;mwrfits, psf, pasta+'\'+'psfgauss.fits'
end