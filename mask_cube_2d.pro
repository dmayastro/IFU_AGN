pro mask_cube_2D

image='C:\Users\incta\Desktop\share\4151_cubes\fits_images\neix.fits'
saida='C:\Users\incta\Desktop\share\4151_cubes\fits_images\neix_M9.fits'

;coordenada x do centro da região circular cujo espectro será extraído
Xc=577.0
;coordenada y do centro da região circular cujo espectro será extraído
Yc=481.0
;raio da região circular cujo espectro será extraído
raio=36.0

extension=0
imagem=MRDFITS(image,extension,header)
z=size(imagem)
imagem_cor=imagem


FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',2
FXADDPAR,headerout,'NAXIS1',z[1]
FXADDPAR,headerout,'NAXIS2',z[2]



   FOR i=0, z[1]-1 DO BEGIN
     FOR j=0, z[2]-1 DO BEGIN
     r=sqrt((i-(Xc-1))^2+(j-(Yc-1))^2)
        IF (r LE raio) THEN BEGIN
           imagem_cor[i,j]=0.0
           print, r
        ENDIF
     ENDFOR
   ENDFOR



MWRFITS, imagem_cor, saida,/CREATE, headerout

END