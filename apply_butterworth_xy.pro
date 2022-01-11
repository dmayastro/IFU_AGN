pro apply_butterworth_xy

pasta = 'F:\AGNs_SINFONI\PSF_correct\psf3368\refrac\'
infile = 'STD_3368_reamost.fits'
filtername = 'STD3368_filter_29_41.fits'

;_____________________________________________________________________________________________________________

cube = mrdfits(pasta+infile, 0, header)
filter = mrdfits(pasta+filtername)
tam = size(cube)
outfile = make_array([tam[1],tam[2],tam[3]],/float,value=0.0)
for k = 0, tam[3] - 1 do outfile[*,*,k] = cube[*,*,k]*filter
header = header[6:*]
mwrfits, outfile, pasta+'fxy_'+infile
end