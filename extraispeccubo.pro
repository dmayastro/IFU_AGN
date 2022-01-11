pro extraispeccubo

;parâmetros de entrada


;caminho completo do cubo de dados partir do qual o espectro da região circular será extraído
entrada='C:\Users\incta\Desktop\m83\2_gas.fits'

;coordenada x do centro da região circular cujo espectro será extraído
;Xc=56.0
;coordenada y do centro da região circular cujo espectro será extraído
;Yc=54.0

;raio da região circular cujo espectro será extraído
;raio=16.0

;arquivo txt que contendo o espectro extraído
saida='C:\Users\incta\Desktop\m83\2_spec.txt'


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

         somaespectros[1,*]=somaespectros[1,*]+cubo[i,j,*]
         ;print, r

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




