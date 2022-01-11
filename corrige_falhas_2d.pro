pro corrige_falhas_2d

image='C:\Users\incta\Desktop\share\4151_cubes\hst\large_oiii_cor26.fits'
saida='C:\Users\incta\Desktop\share\4151_cubes\hst\large_oiii_cor27.fits'

extension=0
imagem=MRDFITS(image,extension,header)
z=size(imagem)
imagem_cor=imagem

;x=32.0

FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',3
FXADDPAR,headerout,'NAXIS1',z[1]
FXADDPAR,headerout,'NAXIS2',z[2]




   ;FOR i=0, z[1]-1 DO BEGIN
    ; FOR j=77, z[2]-1 DO BEGIN
    FOR i=102, 103 DO BEGIN
     FOR j=13, 14 DO BEGIN
       ;IF (imagem[i,j] GE x) THEN BEGIN
           imagem_cor[i,j]=32


     ENDFOR
   ENDFOR





MWRFITS, imagem_cor, saida,/CREATE, headerout

END