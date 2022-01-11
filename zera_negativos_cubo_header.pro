pro zera_negativos_cubo_header

image='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\ESO_sum_qfitsv_e15gig_x25y30n2.fits'
saida='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\ESO_sum_qfitsv_e15gig_x25y30n2_z.fits'

extension=0
imagem=MRDFITS(image,extension,header)
z=size(imagem)
imagem_cor=imagem


lambdazero=FXPAR(header,'CRVAL3')
deltalambda=FXPAR(header,'CDELT3')
reference_pixel=FXPAR(header,'CRPIX3')

FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',3
FXADDPAR,headerout,'NAXIS1',z[1]
FXADDPAR,headerout,'NAXIS2',z[2]
FXADDPAR,headerout,'NAXIS3',z[3]
FXADDPAR,headerout,'CRPIX3',reference_pixel
FXADDPAR,headerout,'CDELT3',deltalambda
FXADDPAR,headerout,'CRVAL3',lambdazero
;FXADDPAR,headerout,'CUNIT3','mum'


FOR k=0, z[3]-1 DO BEGIN
   FOR i=0, z[1]-1 DO BEGIN
     FOR j=0, z[2]-1 DO BEGIN
        IF (imagem[i,j,k] LE 0.0) THEN BEGIN
           imagem_cor[i,j,k]=0.00001
        ENDIF
     ENDFOR
   ENDFOR
ENDFOR




;FOR k=0, z[3]-1 DO BEGIN
;   FOR i=0, z[1]-1 DO BEGIN
;     FOR j=0, z[2]-1 DO BEGIN
;        IF (imagem[i,j,k] LT 0.0) THEN BEGIN
;           imagem_cor[i,j,k]=0.0
;        ENDIF
;     ENDFOR
 ;  ENDFOR
;ENDFOR

FOR k=0, z[3]-1 DO BEGIN
   FOR i=0, z[1]-1 DO BEGIN
     FOR j=0, z[2]-1 DO BEGIN
        IF finite(imagem[i,j,k],/NAN) THEN BEGIN
           imagem_cor[i,j,k]=0.0000001
        ENDIF
     ENDFOR
   ENDFOR
ENDFOR

MWRFITS, imagem_cor, saida,/CREATE, headerout

END