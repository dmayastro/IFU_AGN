pro correct_redshift

;__________________________________________________________________________________________________

entrada = 'C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15.fits'
saida = 'C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15z.fits'
redshift = double(0.01180)
; No caso da velocidade heliocentrica calculada com o rvcorr do IRAF, se vobs=0 ou não existe no header,
; o resultado q ele fornece deve entrar acima com o sinal trocado.
;___________________________________________________________________________________________________

cubo = mrdfits(entrada,0,header)
fxaddpar, header, 'TFIELDS',3
FXBFIND, header, 'CDELT', columns, VALUES_CDELT, N_FOUND
FXBFIND, header, 'CRVAL', columns, VALUES_CRVAL, N_FOUND
;lambdazero = 'CRVAL3'
;deltalambda = 'CDELT3'
lambdazero = 4250.0
deltalambda = 0.691500
tam = size(cubo)
newdelt = make_array(tam[3],/double)
header = header[6:*]

for k = 0.0, tam[3] - 1 do begin
	a = double(lambdazero + deltalambda*k)
	b = double(lambdazero + deltalambda*(k+1))
	lambda1 = double(a/(redshift+1))
	lambda2 = double(b/(redshift+1))
	if k eq 0.0 then begin
		FXADDPAR, HEADER, 'CRVAL3', lambda1
		;astr=STRING(lambda1)
		;stop
		;header[poscrval3] = 'CRVAL3  =             '+astr
	endif
	newdelt[k] = lambda2 - lambda1
endfor
;realnewdelt = double(mean(newdelt))
FXADDPAR, HEADER, 'CDELT3', mean(newdelt)
FXADDPAR, HEADER, 'CD3_3', mean(newdelt)
;realnewdeltstr = string(realnewdelt)
;stop
;header[poscdelt3] = 'CDELT3  =             '+realnewdeltstr
;header[poscd3_3] = 'CD3_3   =             '+realnewdeltstr
;stop
mwrfits, cubo, saida, header
end

