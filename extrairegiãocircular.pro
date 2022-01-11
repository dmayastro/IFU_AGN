pro extrairegiãocircular

;parâmetros de entrada


;caminho completo do cubo de dados partir do qual o espectro da região circular será extraído
entrada='C:\Users\incta\Desktop\share\4151_cubes\images_archives\mask\paper\si_almudena\n4151_velslices\hngc4151_k035_pa0zgrz_x40y40n2g_L.fits'

;coordenada x do centro da região circular cujo espectro será extraído
Xc=60.0
;coordenada y do centro da região circular cujo espectro será extraído
Yc=67.0

;raio da região circular cujo espectro será extraído
raio=3.0

;arquivo txt que contendo o espectro extraído
saida='C:\Users\incta\Desktop\share\4151_cubes\images_archives\mask\paper\si_almudena\n4151_velslices\b1_r3_lower.txt'


;início do programa

cubo=MRDFITS(entrada,0,header)
z=size(cubo)
somaespectros=make_array(2,z[3],/float,value=0.0)

;lambdazero=21045.6
;deltalambda=2.11486
;reference_pixel=1.0
lambdazero=FXPAR(header,'CRVAL3')
deltalambda=FXPAR(header,'CDELT3')
reference_pixel=FXPAR(header,'CRPIX3')

FOR n=1, z[3] DO BEGIN
   somaespectros[0,n-1]=lambdazero+(n-reference_pixel)*deltalambda
ENDFOR

FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      r=sqrt((i-(Xc-1))^2+(j-(Yc-1))^2)
      IF (r LE raio) THEN BEGIN
         somaespectros[1,*]=somaespectros[1,*]+cubo[i,j,*]
         print, r
      ENDIF
   ENDFOR
ENDFOR
somaespectros=somaespectros

OPENW,1,saida
FOR w=0, z[3]-1 DO BEGIN
   PRINTF,1,somaespectros[*,w]
ENDFOR
close,1
print, 'done!'

END




