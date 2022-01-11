pro correct_redshift_v2

;v2 - Colocar somente o redshift da galaxia. Ja corrige automaticamente da velocidade heliocentrica.
;__________________________________________________________________________________________________

pasta = 'C:\Users\Daniel May\Desktop\backup phd\nifs_ngc6951_K_100\
entrada = '127
cubo_header = 'dcsqtedpxrgS20130810S0090


redshift = double(0.00475) ;redshift dentro do double

latitude_obs = 19.823806	;Gemini Norte = 19.823806 e Gemini-Sul = -30.24075
altitude_obs = 4213.		;Gemini-Sul = 2722, Gemini Norte = 4213

;___________________________________________________________________________________________________

cubo = mrdfits(pasta+entrada+'.fits',0,header)
saida = pasta+'vcorr_'+entrada+'.fits'
FXADDPAR, header, 'NAXIS', 3
pix_ref = fxpar(header,'CRPIX3')

deltalambda = fxpar(header,'CD3_3')
lambdazero = fxpar(header,'CRVAL3') + (deltalambda*(1.0 - pix_ref))

tam = size(cubo)
newdelt = make_array(tam[3],/double)
header = header[6:*]
vhel = hel_velocity(pasta+cubo_header+'.fits',latitude_obs,altitude_obs) ;4213m Gemini North Declination=66.10555556

correcao = redshift - (vhel/2.99792458e5)

for k = 0.0, tam[3] - 1 do begin
	a = double(lambdazero + deltalambda*k)
	b = double(lambdazero + deltalambda*(k+1))
	lambda1 = double(a/(correcao+1))
	lambda2 = double(b/(correcao+1))
	if k eq 0.0 then begin
		FXADDPAR, HEADER, 'CRVAL3', lambda1
	endif
	newdelt[k] = lambda2 - lambda1
endfor
;stop
FXADDPAR, HEADER, 'CDELT3', mean(newdelt)
FXADDPAR, HEADER, 'CD3_3', mean(newdelt)
FXADDPAR, HEADER, 'CRPIX3', 1
mwrfits, cubo, saida, header, /create
end

