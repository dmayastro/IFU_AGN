pro Lucy_deconvolution_plus


;input parameters

;name of the data cube (full path) to be deconvolved
cube='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\ESO_sum_qfitsv_e15gig_x25y30n2_add10.fits'

;number of iterations to be used in the deconvolution process
niter=6

;extension of the data cube that will be used
extension=0

;name of the logfile (full path) to be created
logfile='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\logLucy.txt'

;name of the deconvolved data cube (full path)
finalcube='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\ESO_sum_qfitsv_e15gig_x25y30n2_add10_LHalpha.fits'


;beginning of the program

image=MRDFITS(cube,extension,header)
z=size(image)
OPENW,1,logfile
PRINTF,1,'Logfile for lucy_deconvolution_plus'
PRINTF,1,''
PRINTF,1,''
PRINTF,1,'Input file: '+cube
PRINTF,1,'Output file: '+finalcube
PRINTF,1,''

;obtaining the values of lambdazero, deltalambda and reference_pixel
lambdazero=FXPAR(header,'CRVAL3')
deltalambda=FXPAR(header,'CDELT3')
reference_pixel=FXPAR(header,'CRPIX3')

;preparing the header of the deconvolved data cube
FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',3
FXADDPAR,headerout,'NAXIS1',z[1]
FXADDPAR,headerout,'NAXIS2',z[2]
FXADDPAR,headerout,'NAXIS3',z[3]
FXADDPAR,headerout,'CRPIX3',reference_pixel
FXADDPAR,headerout,'CDELT3',deltalambda
FXADDPAR,headerout,'CRVAL3',lambdazero

;enlarging the spatial dimensions
enlarge=''
READ,enlarge,PROMPT='Do you wish to enlarge the spatial dimensions of the datacube before the deconvolution (yes/no)?'
WHILE (enlarge NE 'yes') and (enlarge NE 'no') and (enlarge NE 'y') and (enlarge NE 'n') DO BEGIN
   print, 'ERROR: you must answer yes or no!'
   READ,enlarge,PROMPT='Do you wish to enlarge the spatial dimensions of the datacube before the deconvolution (yes/no)?'
ENDWHILE
IF (enlarge EQ 'yes') or (enlarge EQ 'y') THEN BEGIN
   Nxleft=0
   Nxright=0
   Nydown=0
   Nyup=0
   PRINTF,1,'Spatial dimensions of the datacube enlarged.'
   READ,Nxleft,PROMPT='How many columns do you wish to add at the left side?'
   READ,Nxright,PROMPT='How many columns do you wish to add at the right side?'
   READ,Nydown,PROMPT='How many rows do you wish to add below?'
   READ,Nyup,PROMPT='How many rows do you wish to add above?'
   INxleft=FIX(Nxleft)
   INxright=FIX(Nxright)
   INydown=FIX(Nydown)
   INyup=FIX(Nyup)
   PRINTF,1,INxleft,FORMAT='("Number of spatial pixels added at the left side: ",I)'
   PRINTF,1,INxright,FORMAT='("Number of spatial pixels added at the right side: ",I)'
   PRINTF,1,INydown,FORMAT='("Number of spatial pixels added below: ",I)'
   PRINTF,1,INyup,FORMAT='("Number of spatial pixels added above: ",I)'
   PRINTF,1,''
   newimage=make_array(z[1]+Nxleft+Nxright,z[2]+Nydown+Nyup,z[3],/float)
   FOR i=0, z[1]-1 DO BEGIN
      FOR j=0, z[2]-1 DO BEGIN
         newimage(Nxleft+i,Nydown+j,*)=image(i,j,*)
      ENDFOR
   ENDFOR
   FOR i=0, Nxleft-1 DO BEGIN
      newimage(i,Nydown:Nydown+z[2]-1,*)=image(0,*,*)
   ENDFOR
   FOR i=0, Nxright-1 DO BEGIN
      newimage(Nxleft+z[1]+i,Nydown:Nydown+z[2]-1,*)=image(z[1]-1,*,*)
   ENDFOR
   FOR j=0, Nydown-1 DO BEGIN
      newimage(*,j,*)=newimage(*,Nydown,*)
   ENDFOR
   FOR j=0, Nyup-1 DO BEGIN
      newimage(*,Nydown+z[2]+j,*)=newimage(*,Nydown+z[2]-1,*)
   ENDFOR
ENDIF ELSE BEGIN
   newimage=image
   PRINTF,1,'Spatial dimensions of the datacube were not enlarged.'
   PRINTF,1,''
ENDELSE

;single or variable PSF ?
lucytype=''
READ,lucytype,PROMPT='Do you wish to perform a Richardson-Lucy deconvolution with a single or a variable PSF (single/variable)?'
WHILE (lucytype NE 'single') and (lucytype NE 'variable') and (lucytype NE 's') and (lucytype NE 'v') DO BEGIN
   print, 'ERROR: the type of the deconvolution must be single or variable!'
   READ,lucytype,PROMPT='Do you wish to perform a Richardson-Lucy deconvolution with a single or a variable PSF (single/variable)?'
ENDWHILE

;deconvolutions with a single PSF
IF (lucytype EQ 'single') or (lucytype EQ 's') THEN BEGIN
   PSFtype=''
   READ,PSFtype,PROMPT='Do you want to use a gaussian or a real PSF (gaussian/real)?'
   WHILE (PSFtype NE 'gaussian') and (PSFtype NE 'real') and (PSFtype NE 'g') and (PSFtype NE 'r') DO BEGIN
      print, 'ERROR: you must choose a real or a gaussian PSF!'
      READ,PSFtype,PROMPT='Do you want to use a gaussian or a real PSF (gaussian/real)?
   ENDWHILE

   ;deconvolution with a gaussian and single PSF
   IF (PSFtype EQ 'gaussian') or (PSFtype EQ 'g') THEN BEGIN
      PRINTF,1,'The process was performed with a single and gaussian PSF.'
      READ,FWHM,PROMPT='Enter the value of the FWHM (in pixels) to be used to construct a gaussian PSF image:'
      PRINTF,1,FWHM,FORMAT='("Value of the FWHM (in pixels) used to create the gaussian PSF: ",F)'
      deconv=newimage
      z=size(newimage)
      psf=psf_Gaussian(NPIXEL=[z[1],z[2]], FWHM=FWHM, /NORMALIZE)
      FOR n=0, z[3]-1 DO BEGIN
         FOR k=1, niter do Max_Likelihood, newimage(*,*,n), psf, deconvlucy
         deconv(*,*,n)=deconvlucy
         deconvlucy=0
         print, n+1
      ENDFOR
   ENDIF

   ;deconvolution with a real and single PSF
   IF (PSFtype EQ 'real') or (PSFtype EQ 'r') THEN BEGIN
      PRINTF,1,'The process was performed with a single and real PSF.'
      PSF=''
      READ,PSF,PROMPT='Enter the name (full path) of the image to be used as PSF in the deconvolution:'
      PRINTF,1,'File used as PSF: '+PSF
      psfsemnorm=MRDFITS(PSF)
      deconv=newimage
      size=size(newimage)
      soma=0.0
      z=size(psfsemnorm)
      FOR i=0, z[1]-1 DO BEGIN
         FOR j=0, z[2]-1 DO BEGIN
            soma=soma+psfsemnorm(i,j)
         ENDFOR
      ENDFOR
      psf=psfsemnorm/soma
      FOR n=0, size[3]-1 DO BEGIN
         FOR k=1, niter do Max_Likelihood, newimage(*,*,n), psf, deconvlucy
         deconv(*,*,n)=deconvlucy
         deconvlucy=0
         print, n+1
      ENDFOR
   ENDIF
ENDIF

;deconvolutions with a variable PSF
IF (lucytype EQ 'variable') or (lucytype EQ 'v') THEN BEGIN
   PRINTF,1,'The process was performed with a variable PSF.'
   method=''
   READ,method,PROMPT='Which method do you want to use to determine the PSF for the Richardson-Lucy deconcolution (1/2/3)?'
   WHILE (method NE '1') and (method NE '2') and (method NE '3') DO BEGIN
      print, 'ERROR: you must choose between methods 1, 2 and 3!'
      READ,method,PROMPT='Which method do you want to use to determine the PSF for the Richardson-Lucy deconcolution (1/2/3)?'
   ENDWHILE

   ;deconvolution with method 1
   IF (method EQ '1') THEN BEGIN
      unity=''
      PRINTF,1,'Method used to determine the PSF for the deconvolution: '+method
      READ,sizepix,PROMPT='Enter the size of each spatial pixel in arcsec (sizepix):'
      READ,lambdaref,PROMPT='Enter the value of the reference wavelength (lambdaref):'
      READ,unity,PROMPT='Enter the unity of the reference wavelength:'
      READ,FWHMref,PROMPT='Enter the value of the reference FWHM (FWHMref) in pixels:'
      PRINTF,1,sizepix,FORMAT='("Size of each spatial pixel in arcsec (sizepix): ",F)'
      PRINTF,1,lambdaref,'("Value of the reference wavelength (lambdaref): ",F)'
      PRINTF,1,'Unity used for wavelength measure: '+unity
      PRINTF,1,FWHMref,FORMAT='("Value of the reference FWHM (FWHMref) in pixels: ",F)'
      z = size(newimage)
      deconv = MAKE_ARRAY(z[1], z[2], z[3], /float)
      FOR i = 0, z[3] - 1 DO BEGIN
         n=i+1
         lambda = lambdazero+(n-reference_pixel)*deltalambda
         lma = (FWHMref*sizepix)*(lambda^(-0.484))/(lambdaref^(-0.484))
         lma = lma / sizepix
         psf = psf_gaussian(NPIXEL = [z[1],z[2]], FWHM = lma, $
         /normalize)
         FOR k = 1, niter DO Max_Likelihood, newimage[*,*,i], $
         psf, deconvlucy
         deconv[*,*,i] = deconvlucy
         ;psf = 0
         deconvlucy = 0
         print , n
      ENDFOR
   ENDIF

   ;deconvolution with method 2
   IF (method EQ '2') THEN BEGIN
      unity=''
      PRINTF,1,'Method used to determine the PSF for the deconvolution: '+method
      READ,FWHMB,PROMPT='Enter the value of the FWHM at the blue border of the spectrum (FWHMB) in pixels:'
      READ,FWHMR,PROMPT='Enter the value of the FWHM at the red border of the spectrum (FWHMR) in pixels:'
      lambdaB=lambdazero+(1-reference_pixel)*deltalambda
      lambdaR=lambdazero+(z[3]-reference_pixel)*deltalambda
      READ,unity,PROMPT='Enter the unity of wavelengths provided above:'
      PRINTF,1,FWHMB,FORMAT='("Value of the FWHM at the blue border of the spectrum (FWHMB) in pixels: ",F)'
      PRINTF,1,FWHMR,FORMAT='("Value of the FWHM at the red border of the spectrum (FWHMR) in pixels: ",F)'
      PRINTF,1,lambdaB,FORMAT='("Value of the wavelength at the blue border of the spectrum (lambdaB): ",F)'
      PRINTF,1,lambdaR,FORMAT='("Value of the wavelength at the red border of the spectrum (lambdaR): ",F)'
      PRINTF,1,'Unity used for wavelength measure: '+unity
      deconv=newimage
      z=size(newimage)
      psfB=psf_Gaussian(NPIXEL=[z[1],z[2]], FWHM=FWHMB, /NORMALIZE)
      psfR=psf_Gaussian(NPIXEL=[z[1],z[2]], FWHM=FWHMR, /NORMALIZE)
      FOR y=0, z[3]-1 DO BEGIN
         n=y+1
         lambda=lambdazero+(n-reference_pixel)*deltalambda
         psf=psfB*(1-((lambda-lambdaB)/(lambdaR-lambdaB)))+psfR*((lambda-lambdaB)/(lambdaR-lambdaB))
         FOR k=1, niter do Max_Likelihood, newimage(*,*,y), psf, deconvlucy
         deconv(*,*,y)=deconvlucy
         deconvlucy=0
         print, n
      ENDFOR
   ENDIF

   ;deconvolution with method 3
   IF (method EQ '3') THEN BEGIN
      unity=''
      PRINTF,1,'Method used to determine the PSF for the deconvolution: '+method
      PSFB=''
      PSFR=''
      READ,PSFB,PROMPT='Enter the name (full path) of the image representing the PSF at the blue border of the spectrum (PSFB):'
      READ,PSFR,PROMPT='Enter the name (full path) of the image representing the PSF at the red border of the spectrum (PSFR):'
      lambdaB=lambdazero+(1-reference_pixel)*deltalambda
      lambdaR=lambdazero+(z[3]-reference_pixel)*deltalambda
      READ,unity,PROMPT='Enter the unity of wavelengths provided above:'
      PRINTF,1,'Name (full path) of the image representing the PSF at the blue border of the spectrum (PSFB): '+PSFB
      PRINTF,1,'Name (full path) of the image representing the PSF at the red border of the spectrum (PSFR): '+PSFR
      PRINTF,1,lambdaB,FORMAT='("Value of the wavelength at the blue border of the spectrum (lambdaB): ",F)'
      PRINTF,1,lambdaR,FORMAT='("Value of the wavelength at the red border of the spectrum (lambdaR): ",F)'
      PRINTF,1,'Unity used for wavelength measure: '+unity
      psfblue_semnorm=MRDFITS(PSFB)
      psfred_semnorm=MRDFITS(PSFR)
      z=size(psfblue_semnorm)
      soma=0.0
      FOR i=0, z[1]-1 DO BEGIN
         FOR j=0, z[2]-1 DO BEGIN
            soma=soma+psfblue_semnorm(i,j)
         ENDFOR
      ENDFOR
      psfblue=psfblue_semnorm/soma
      soma=0.0
      FOR i=0, z[1]-1 DO BEGIN
         FOR j=0, z[2]-1 DO BEGIN
            soma=soma+psfred_semnorm(i,j)
         ENDFOR
      ENDFOR
      psfred=psfred_semnorm/soma
      deconv=newimage
      size=size(newimage)
      FOR y=0, size[3]-1 DO BEGIN
         n=y+1
         lambda=lambdazero+(n-reference_pixel)*deltalambda
         psf=psfblue*(1-((lambda-lambdaB)/(lambdaR-lambdaB)))+psfred*((lambda-lambdaB)/(lambdaR-lambdaB))
         FOR k=1, niter do Max_Likelihood, newimage(*,*,y), psf, deconvlucy
         deconv(*,*,y)=deconvlucy
         deconvlucy=0
         print, n
      ENDFOR
   ENDIF
ENDIF

IF (enlarge EQ 'yes') or (enlarge EQ 'y') THEN BEGIN
   z=size(deconv)
   final_image=deconv(Nxleft:z[1]-1-Nxright,Nydown:z[2]-1-Nyup,*)
ENDIF ELSE BEGIN
   final_image=deconv
ENDELSE
close,1
MWRFITS, final_image, finalcube, headerout
print, "Done!"


END
