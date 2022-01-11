;Cria um cubo artificial qualquer
pro cube_create

;muda o centro da gaussiana no espectro tam_z, em pixeis
mi = -5.0

tam_z = 100 ;fixo para um cubo, e para todos os outros q se somarem a ele.

fwhm=7.0
;muda o centro do blob
xc = 20
yc = 20

raio = 10


saida = '\\DMAY-PC\Users\dmay\Desktop\sintetic_cube_fwhm_7_mi_m5.fits'
z=make_array(3,/INT)
z[0]=60
z[1]=60
z[2]=tam_z
cubo=make_array(z[0],z[1],z[2],/FLOAT)
extension=0

lambdazero=21168.0
deltalambda=2.0
reference_pixel=1.0

FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',3
FXADDPAR,headerout,'NAXIS1',z[0]
FXADDPAR,headerout,'NAXIS2',z[1]
FXADDPAR,headerout,'NAXIS3',z[2]
FXADDPAR,headerout,'CRPIX3',reference_pixel
FXADDPAR,headerout,'CDELT3',deltalambda
FXADDPAR,headerout,'CRVAL3',lambdazero

ln2=0.693147  ;0.693147180
sigma = fwhm / (2*sqrt(2*ln2))
pi=3.14159  ;3.1415926535
k1= 1/sqrt(2*pi*sigma^2)
k2= 2*(sigma^2)

;espectro = make_array(tam_z, /double)
;espectro normalizado:
;if tipo eq '1' then begin
  FOR i=0, z[0]-1 DO BEGIN
     FOR j=0, z[1]-1 DO BEGIN
       for k=0, tam_z-1 do begin
       v1 = (-1)*(k-(tam_z/2)-mi)^2
       valor = k1 * exp( v1/k2 )
      ;valor = exp( v1/k2 ) ;para espectro não normalizado
         r=sqrt((i-(xc-1))^2+(j-(yc-1))^2)
         IF (r LE raio) THEN BEGIN
           cubo[i,j,k] = valor*(raio+1-r)^2 ;cria blob
           ;print, r
         ENDIF ELSE BEGIN
           cubo[i,j,k] = valor
	     ENDELSE
	   endfor
     ENDFOR
  ENDFOR
;endif

mwrfits,cubo,saida,headerout,/CREATE

print,"Done!"

close,/ALL

end
