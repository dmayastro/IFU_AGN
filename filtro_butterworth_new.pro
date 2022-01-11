pro filtro_Butterworth_new

;parâmetros de entrada

;nome do cubo de dados
cube='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zcgi.fits'

;extensão do cubo de dados que será utilizado no processo
extension=0

cutoff_frequency_x=0.25

cutoff_frequency_y=0.30

order=2.0

elipse='no'
elipse2='no'
elipsevssquare='yes'

;nome do cubo final
saida='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zcgi_x25y30n2.fits'

;output FFT image (padded, image taken in the middle of z-axis)
namefft1="C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\fft.fits"
;output FFT*filter image (padded, image taken in the middle of z-axis))
namefft2="C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\fft_filtered.fits"
;########################################################################


;início do programa

image=MRDFITS(cube,extension,header)
z=size(image)
newimage=make_array(z[1]+6,z[2]+6,z[3],/float)
FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      newimage(3+i,3+j,*)=image(i,j,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      newimage(i,3+j,*)=image(0,j,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      newimage(3+z[1]+i,3+j,*)=image(z[1]-1,j,*)
   ENDFOR
ENDFOR
FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(3+i,j,*)=image(i,0,*)
   ENDFOR
ENDFOR
FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(3+i,3+z[2]+j,*)=image(i,z[2]-1,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(i,j,*)=newimage(3,j,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(3+z[1]+i,j,*)=newimage(3+z[1]-1,j,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(i,3+z[2]+j,*)=newimage(3,3+z[2]+j,*)
   ENDFOR
ENDFOR
FOR i=0, 3-1 DO BEGIN
   FOR j=0, 3-1 DO BEGIN
      newimage(3+z[1]+i,3+z[2]+j,*)=newimage(3+z[1]-1,3+z[2]+j,*)
   ENDFOR
ENDFOR

z=size(newimage)
cubo_final=make_array(z[1]+2*15,z[2]+2*15,z[3],/FLOAT)

cubo_final[15:z[1]+15-1,15:z[2]+15-1,*]=newimage
interpol_values=make_array(2,2,/FLOAT)
temporary_column=make_array(2,z[2]+2,/FLOAT)
temporary_column[0,0]=0.0
temporary_column[0,z[2]+1]=z[2]+2*15-1
temporary_column[0,1:z[2]]=15+FINDGEN(z[2])
temporary_column[1,0]=0.0
temporary_column[1,z[2]+1]=0.0
FOR k=0, z[3]-1 DO BEGIN
   interpol_values[1,0]=0.0
   interpol_values[0,0]=0.0
   interpol_values[0,1]=15
   FOR i=15, z[1]+15-1 DO BEGIN
      interpol_values[1,1]=cubo_final[i,15,k]
      cubo_final[i,0:15,k]=INTERPOL(interpol_values[1,*], interpol_values[0,*], FINDGEN(15+1))
   ENDFOR

   interpol_values[0,0]=z[2]+15-1
   interpol_values[0,1]=z[2]+2*15-1
   interpol_values[1,1]=0.0
   FOR i=15, z[1]+15-1 DO BEGIN
      interpol_values[1,0]=cubo_final[i,z[2]+15-1,k]
      cubo_final[i,z[2]+15-1:z[2]+2*15-1,k]=INTERPOL(interpol_values[1,*], interpol_values[0,*], z[2]+15-1+FINDGEN(15+1))
   ENDFOR

   interpol_values[0,0]=0.0
   interpol_values[0,1]=15
   interpol_values[1,0]=0.0
   FOR j=15, z[2]+15-1 DO BEGIN
      interpol_values[1,1]=cubo_final[15,j,k]
      cubo_final[0:15,j,k]=INTERPOL(interpol_values[1,*], interpol_values[0,*], FINDGEN(15+1))
   ENDFOR

   interpol_values[0,0]=z[1]+15-1
   interpol_values[0,1]=z[1]+2*15-1
   interpol_values[1,1]=0.0
   FOR j=15, z[2]+15-1 DO BEGIN
      interpol_values[1,0]=cubo_final[z[1]+15-1,j,k]
      cubo_final[z[1]+15-1:z[1]+2*15-1,j,k]=INTERPOL(interpol_values[1,*], interpol_values[0,*], z[1]+15-1+FINDGEN(15+1))
   ENDFOR

   FOR i=0, 15-1 DO BEGIN
      temporary_column[1,1:z[2]]=cubo_final[i,15:z[2]+15-1,k]
      cubo_final[i,*,k]=INTERPOL(temporary_column[1,*], temporary_column[0,*], FINDGEN(z[2]+2*15))
   ENDFOR

   FOR i=z[1]+15, z[1]+2*15-1 DO BEGIN
      temporary_column[1,1:z[2]]=cubo_final[i,15:z[2]+15-1,k]
      cubo_final[i,*,k]=INTERPOL(temporary_column[1,*], temporary_column[0,*], FINDGEN(z[2]+2*15))
   ENDFOR

   FOR j=0, 15-1 DO BEGIN
      cubo_final[*,j,k]=INTERPOL(cubo_final[*,j,k], FINDGEN(z[1]+2*15), FINDGEN(z[1]+2*15))
   ENDFOR

   FOR j=z[2]+15, z[2]+2*15-1 DO BEGIN
      cubo_final[*,j,k]=INTERPOL(cubo_final[*,j,k], FINDGEN(z[1]+2*15), FINDGEN(z[1]+2*15))
   ENDFOR

   print, k

ENDFOR

z=SIZE(cubo_final)

imfilter=MAKE_ARRAY(2*z[1],2*z[2],/DOUBLE)

aa=cutoff_frequency_x*z[1]
bb=cutoff_frequency_y*z[2]
cc=cutoff_frequency_x*z[1]
dd=cutoff_frequency_y*z[2]

exporder=2.0*order
;constructs of the filter image 1
imfilter1=imfilter
;xm=z[1]-0.5
;ym=z[2]-0.5
xm=z[1]
ym=z[2]


FOR i=0, 2*z[1]-1 DO BEGIN
   FOR j=0, 2*z[2]-1 DO BEGIN
      duv=SQRT(((i-xm)/aa)^2.0+((j-ym)/bb)^2.0)
      argfrac=duv^exporder
      imfilter1(i,j)=1.0/(1.0+argfrac)
   ENDFOR
ENDFOR
;constructs of the filter image 2
imfilter2=imfilter
FOR i=0, 2*z[1]-1 DO BEGIN
   FOR j=0, 2*z[2]-1 DO BEGIN
      duv1=ABS(i-xm)/cc
      duv2=ABS(j-ym)/dd
      argfrac1=duv1^exporder
      argfrac2=duv2^exporder
      aux1=1.0/(1.0+argfrac1)
      aux2=1.0/(1.0+argfrac2)
      imfilter2(i,j)=aux1*aux2
   ENDFOR
ENDFOR

IF (elipse EQ 'yes') THEN BEGIN
   imfilter=imfilter1
ENDIF

IF (elipse2 EQ 'yes') THEN BEGIN
   imfilter=imfilter1*imfilter1
ENDIF

IF (elipsevssquare EQ 'yes') THEN BEGIN
   imfilter=imfilter1*imfilter2
ENDIF

cubeout=cubo_final
npadd=MAKE_ARRAY(2*z[1],2*z[2],/DOUBLE)
imfft1=npadd
imfft2=npadd

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

k=0
zc=FLOOR(z[3]/2.0)
;starts the filtering
FOR k=0, z[3]-1 DO BEGIN
   imagein=cubo_final(*,*,k)
   imageflag=npadd
   ;makes the 'padding' of each image
   FOR i=0, 2*z[1]-1 DO BEGIN
      FOR j=0, 2*z[2]-1 DO BEGIN
         IF ((i LE z[1]-1) AND (j GT z[2]-1)) THEN BEGIN
            flag=((i+j) MOD 2) NE 0 ? 1 : -1
            imageflag(i,j)=flag*imagein(i,j-z[2])
         ENDIF ELSE BEGIN
            imageflag(i,j)=0
         ENDELSE
      ENDFOR
   ENDFOR

   twodfft=FFT(imageflag)
   IF (k EQ zc) THEN imfft1=ABS(twodfft) ;FFT image (modulus), just for check
   guv=twodfft*imfilter
   IF (k EQ zc) THEN imfft2=ABS(guv)  ;FFT*FILTER image (modulus), just for check
   invfft=REAL_PART(FFT(guv,/INVERSE))

   paddout=npadd
   imageout=imagein
   ;crops the final image
   FOR i=0, 2*z[1]-1 DO BEGIN
      FOR j=0, 2*z[2]-1 DO BEGIN
         flag=((i+j) MOD 2) NE 0 ? 1 : -1
         paddout(i,j)=flag*invfft(i,j)
         IF ((i LE z[1]-1) AND (j GT z[2]-1)) THEN imageout(i,(j-z[2]))=paddout(i,j)
      ENDFOR
   ENDFOR
   cubeout(*,*,k)=imageout
   print, k+1
ENDFOR

cubeout_final=cubeout[18:z[1]-19,18:z[2]-19,*]

MWRFITS,cubeout_final,saida,headerout,/CREATE
MWRFITS,imfft1,namefft1,headerout,/CREATE
MWRFITS,imfft2,namefft2,headerout,/CREATE
print,"Done!"

END





