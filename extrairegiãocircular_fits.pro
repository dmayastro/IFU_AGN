pro extrairegiãocircular_fits

;parâmetros de entrada


;caminho completo do cubo de dados partir do qual o espectro da região circular será extraído
entrada='C:\Users\incta\Desktop\share\4151_cubes\4151_banda-k-sivii\kl_si2_i_x35y40n2_L10T1z_ESSE.fits'

;coordenada x do centro da região circular cujo espectro será extraído
Xc=73.0
;coordenada y do centro da região circular cujo espectro será extraído
Yc=66.0

;raio da região circular cujo espectro será extraído
raio=5.0

;arquivo txt que contendo o espectro extraído
saida='C:\Users\incta\Desktop\share\4151_cubes\4151_banda-k-sivii\blob_west_r6_si7_b.fits'


;início do programa

cubo=MRDFITS(entrada,0,header)
z=size(cubo)
somaespectros=make_array(z[3],/float,value=0.0)

lambdazero=FXPAR(header,'CRVAL3')
deltalambda=FXPAR(header,'CDELT3')
reference_pixel=FXPAR(header,'CRPIX3')

;lambdazero=11471.1
;deltalambda=1.05179
;reference_pixel=1.0

;preparing the header of the deconvolved data cube
FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS', 1
FXADDPAR,headerout,'DISPAXIS',1
FXADDPAR,headerout,'CTYPE1','LINEAR'
FXADDPAR,headerout,'NAXIS1',z[3]
FXADDPAR,headerout,'CRPIX1',reference_pixel
FXADDPAR,headerout,'CDELT1',deltalambda
FXADDPAR,headerout,'CRVAL1',lambdazero
FXADDPAR,headerout,'WAT1_001','wtype=linear label=Wavelength units=angstroms units_display=angstrom'
FXADDPAR,headerout,'WAT1_002','s'
FXADDPAR,headerout, 'CUNIT1', 'angstroms'
FXADDPAR,headerout, 'CUNIT2', 'angstroms'
FXADDPAR,headerout, 'CUNIT3', 'angstroms'



;FOR n=1, z[3] DO BEGIN
;   somaespectros[0,n-1]=lambdazero+(n-reference_pixel)*deltalambda
;ENDFOR

FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      r=sqrt((i-(Xc-1))^2+(j-(Yc-1))^2)
      IF (r LE raio) THEN BEGIN
         somaespectros=somaespectros+cubo[i,j,*]
         print, r
      ENDIF
   ENDFOR
ENDFOR
somaespectros=somaespectros
;______________________________quadrado-todo cubo:

;_____________________________

;OPENW,1,saida
;FOR w=0, z[3]-1 DO BEGIN
;   PRINTF,1,somaespectros[*,w]
;ENDFOR
;close,1

MWRFITS,somaespectros,saida,headerout,/CREATE


print, 'done!'

END