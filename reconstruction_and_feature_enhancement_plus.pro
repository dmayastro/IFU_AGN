pro reconstruction_and_feature_enhancement_plus


;input parameters

;directory containing ALL files that will be used during the execution of the program
infolder='C:\Users\incta\Desktop\share\Hektor_VIMOS\CALIBRATED-abell14\new\pca_hbo_red_linemask\'

;name (full path) of the reconstructed data cube (full path)
finalcube='C:\Users\incta\Desktop\share\Hektor_VIMOS\CALIBRATED-abell14\new\pca_hbo_red_linemask\no234.fits'

;name of the logfile (full path) to be created
logfile='C:\Users\incta\Desktop\share\Hektor_VIMOS\CALIBRATED-abell14\new\pca_hbr_blue_linemask\log_reconstructed.txt'

;values of the lowest and highest spectral pixel used in the execution of the PCA. If the user has made a re-binning
;in the data cube before the execution of the PCA, then, the parameters minspecpx and maxspecpx must correspond to
;the values of the already binned data cube
minspecpx=1
maxspecpx=2580

;values of the lowest and highest spatial pixel of the horizontal axis used in the execution of the PCA
minspatpx_x=1
maxspatpx_x=66

;values of the lowest and highest spatial pixel of the vertical axis used in the execution of the PCA
minspatpx_y=1
maxspatpx_y=64

;value of the wavelength corresponding to the lowest spectral pixel (minspecpx) considered in the execution of the
;PCA Tomography
lambdazero=5701.4

;wavelength interval between two consecutive spectral pixels
deltalambda=0.62


;beginning of the program

m=maxspecpx-minspecpx+1
n=(maxspatpx_x-minspatpx_x+1)*(maxspatpx_y-minspatpx_y+1)
final_cube=make_array(maxspatpx_x-minspatpx_x+1,maxspatpx_y-minspatpx_y+1,m,/FLOAT)
Ilinha_beta_lambda=make_array(m,n,/FLOAT)
method=''
reconstruction_type=''
OPENW,2,logfile
PRINTF,2,'Logfile for feature_enhancement_plus'
PRINTF,2,''
PRINTF,2,''
PRINTF,2,'Output file: '+finalcube
PRINTF,2,''

;preparing the header of the corrected data cube
FXADDPAR,headerout,'SIMPLE','T'
FXADDPAR,headerout,'BITPIX',-32
FXADDPAR,headerout,'NAXIS',3
FXADDPAR,headerout,'NAXIS1',maxspatpx_x-minspatpx_x+1
FXADDPAR,headerout,'NAXIS2',maxspatpx_y-minspatpx_y+1
FXADDPAR,headerout,'NAXIS3',m
FXADDPAR,headerout,'CRPIX3',981
FXADDPAR,headerout,'CDELT3',deltalambda
FXADDPAR,headerout,'CRVAL3',lambdazero


;simple reconstruction or feature suppression and enhancement ?
READ,reconstruction_type,PROMPT='Do you wish to perform a simple reconstruction of the datacube or a feature suppression and enhancement process (reconstruction/feature)?'
WHILE (reconstruction_type NE 'reconstruction') and (reconstruction_type NE 'r') and (reconstruction_type NE 'feature') and (reconstruction_type NE 'f') DO BEGIN
   print, 'ERROR: you must answer "reconstruction" or "feature"!'
   READ,reconstruction_type,PROMPT='Do you wish to perform a simple reconstruction of the datacube or a feature suppression and enhancement process (reconstruction/feature)?'
ENDWHILE

;simple reconstruction process
IF (reconstruction_type EQ 'reconstruction') or (reconstruction_type EQ 'r') THEN BEGIN
   PRINTF,2,'A simple reconstruction process was applied.'
   PRINTF,2,''
   READ,method,PROMPT='Do you wish to execute the process using the constructed eigenspectra and tomograms or using the tables containing the parameters of the eigenspectra and tomograms (eigenspectra/tables)?'
   WHILE (method NE 'eigenspectra') and (method NE 'e') and (method NE 'tables') and (method NE 't') DO BEGIN
      print, 'ERROR: you must answer "eigenspectra" or "tables"!'
      READ,method,PROMPT='Do you wish to execute the process using the constructed eigenspectra and tomograms or using the tables containing the parameters of the eigenspectra and tomograms (eigenspectra/tables)?'
   ENDWHILE

   ;reconstruction using the eigenspectra and tomograms
   IF (method EQ 'eigenspectra') or (method EQ 'e') THEN BEGIN
      eigenspectra_prefix=''
      tomograms_prefix=''
      r=0
      READ,eigenspectra_prefix,PROMPT='Enter the prefix of the eigenspectra files: '
      READ,tomograms_prefix,PROMPT='Enter the prefix of the tomograms files: '
      READ,r,PROMPT='Enter the number of the highest order eigenvector to be used in the reconstruction process: '
      T_beta_k=make_array(r,n,/FLOAT)
      E_lambda_k=make_array(r,m,/FLOAT)
      eigenspectrum=make_array(2,m,/FLOAT)
      tomogram=make_array(maxspatpx_x-minspatpx_x+1,maxspatpx_y-minspatpx_y+1,/FLOAT)
      FOR w=0, r-1 DO BEGIN
         l=STRING(w+1)
         l=STRTRIM(l,1)
         OPENR,1,infolder+eigenspectra_prefix+l+'.txt'
         READF,1,eigenspectrum
         close,1
         E_lambda_k[w,*]=eigenspectrum[1,*]
         tomogram=MRDFITS(infolder+tomograms_prefix+l+'.fits')
         y=0
         FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
            FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
               T_beta_k[w,y]=tomogram[i,j]
               y=y+1
            ENDFOR
         ENDFOR
      ENDFOR
      Ilinha_beta_lambda=T_beta_k##TRANSPOSE(E_lambda_k)
   ENDIF

   ;reconstruction using the tables
   IF (method EQ 'tables') or (method EQ 't') THEN BEGIN
      SCORE=''
      PC=''
      r=0
      T_beta_k=make_array(m,n,/FLOAT)
      E_lambda_k=make_array(m,m,/FLOAT)
      READ,SCORE,PROMPT="Enter the file's name of the table containing the parameters of the tomograms obtained with the PCA tomography: "
      READ,PC,PROMPT="Enter the file's name of the table containing the parameters of the eigenvectors obtained with the PCA tomography: "
      OPENR,1,infolder+SCORE
      READF,1,T_beta_k
      close,1
      OPENR,1,infolder+PC
      READF,1,E_lambda_k
      close,1
      READ,r,PROMPT='Enter the number of the highest order eigenvector to be used in the reconstruction process: '
      FOR w=r, m-1 DO BEGIN
         E_lambda_k[w,*]=E_lambda_k[w,*]*0.0
      ENDFOR
      Ilinha_beta_lambda=T_beta_k##TRANSPOSE(E_lambda_k)
   ENDIF
   rstring=STRING(r)
   rstring=STRTRIM(r,1)
   PRINTF,2,'The first '+rstring+' eigenvectors were used in the reconstruction process.'
   PRINTF,2,''

   ;addition of the average spectrum
   mean=''
   READ,mean,PROMPT='Do you wish to add the average spectrum to the reconstructed datacube (yes/no)?'
   WHILE (mean NE 'yes') and (mean NE 'y') and (mean NE 'no') and (mean NE 'n') DO BEGIN
      print, 'ERROR: you must answer "yes" or "no"!'
      READ,mean,PROMPT='Do you wish to add the average spectrum to the reconstructed datacube (yes/no)?'
   ENDWHILE
   IF (mean EQ 'yes') or (mean EQ 'y') THEN BEGIN
      PRINTF,2,'The average spectrum was added to the final datacube.'
      Q_lambda=make_array(2,m,/FLOAT)
      average_spectrum=''
      Ilinha_beta_lambda_zero=make_array(m,n,/FLOAT)
      READ,average_spectrum,PROMPT='Enter the name of the file containing the average spectrum of the original datacube used in the PCA tomography: '
      OPENR,1,infolder+average_spectrum
      READF,1,Q_lambda
      close,1
      FOR lambda=0, m-1 DO BEGIN
         FOR beta=0, n-1 DO BEGIN
            Ilinha_beta_lambda_zero[lambda,beta]=Ilinha_beta_lambda[lambda,beta]+Q_lambda[1,lambda]
         ENDFOR
      ENDFOR
      y=0
      FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
         FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
            final_cube[i,j,*]=Ilinha_beta_lambda_zero[*,y]
            y=y+1
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      PRINTF,2,'The average spectrum was not added to the final datacube.'
      y=0
      FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
         FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
            final_cube[i,j,*]=Ilinha_beta_lambda[*,y]
            y=y+1
         ENDFOR
      ENDFOR
   ENDELSE
   MWRFITS,final_cube,finalcube,headerout
ENDIF

;feature suppression and enhancement process
IF (reconstruction_type EQ 'feature') or (reconstruction_type EQ 'f') THEN BEGIN
   PRINTF,2,'A feature suppresion and enhancement proccess was applied.'
   PRINTF,2,''
   READ,method,PROMPT='Do you wish to execute the process using the constructed eigenspectra and tomograms or using the tables containing the parameters of the eigenspectra and tomograms (eigenspectra/tables)?'
   WHILE (method NE 'eigenspectra') and (method NE 'e') and (method NE 'tables') and (method NE 't') DO BEGIN
      print, 'ERROR: you must answer "eigenspectra" or "tables"!'
      READ,method,PROMPT='Do you wish to execute the process using the constructed eigenspectra and tomograms or using the tables containing the parameters of the eigenspectra and tomograms (eigenspectra/tables)?'
   ENDWHILE

   ;feature suppression and enhancement executed using the eigenspectra and tomograms
   IF (method EQ 'eigenspectra') or (method EQ 'e') THEN BEGIN
      eigenspectra_prefix=''
      tomograms_prefix=''
      feature_method=''
      READ,eigenspectra_prefix,PROMPT='Enter the prefix of the eigenspectra files: '
      READ,tomograms_prefix,PROMPT='Enter the prefix of the tomograms files: '
      READ,feature_method,PROMPT='Do you wish to perform the feature suppression and enhancement process using method 1 or method 2 (1/2)? '
      WHILE (feature_method NE '1') and (feature_method NE '2') DO BEGIN
         print, 'ERROR: you must answer "1" or "2"!'
         READ,feature_method,PROMPT='Do you wish to perform the feature suppression and enhancement process using method 1 or method 2 (1/2)?'
      ENDWHILE

      ;method 1 for the feature suppression and enhancement process
      IF (feature_method EQ '1') THEN BEGIN
         PRINTF,2,'Method used for the feature suppression and enhancement procedure: '+feature_method
         PRINTF,2,''
         gama_k_list=''
         eigenvectors_list=''
         READ,eigenvectors_number,PROMPT='Enter the number of the eigenvectors to be considered in the feature suppression and enhancement process: '
         READ,eigenvectors_list,PROMPT='Enter the name of the file containing the eigenvectors to be considered in the process: '
         READ,gama_k_list,PROMPT='Enter the name of the file containing the feature factors to be attributed to all the eigenvectors in the process: '
         gama_k=make_array(eigenvectors_number,/FLOAT)
         eigenvectors=make_array(eigenvectors_number,/INTEGER)
         OPENR,1,infolder+gama_k_list
         READF,1,gama_k
         close,1
         OPENR,1,infolder+eigenvectors_list
         READF,1,eigenvectors
         close,1
         PRINTF,2,'Eigenvectors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,eigenvectors[w]
         ENDFOR
         PRINTF,2,''
         PRINTF,2,'Feature factors associated to each one of the above eigenvectors: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,gama_k[w]
         ENDFOR
         PRINTF,2,''
         E_lambda_k=make_array(eigenvectors_number,m,/FLOAT)
         T_beta_k=make_array(eigenvectors_number,n,/FLOAT)
         E_lambda_k_gama=make_array(eigenvectors_number,m,/FLOAT)
         eigenspectrum=make_array(2,m,/FLOAT)
         tomogram=make_array(maxspatpx_x-minspatpx_x+1,maxspatpx_y-minspatpx_y+1,/FLOAT)
         FOR w=0, eigenvectors_number-1 DO BEGIN
            l=STRING(eigenvectors[w])
            l=STRTRIM(l,1)
            OPENR,1,infolder+eigenspectra_prefix+l+'.txt'
            READF,1,eigenspectrum
            close,1
            E_lambda_k[w,*]=eigenspectrum[1,*]
            tomogram=MRDFITS(infolder+tomograms_prefix+l+'.fits')
            y=0
            FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
               FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
                  T_beta_k[w,y]=tomogram[i,j]
                  y=y+1
               ENDFOR
            ENDFOR
         ENDFOR
         FOR w=0, eigenvectors_number-1 DO BEGIN
            E_lambda_k_gama[w,*]=E_lambda_k[w,*]*gama_k[w]
         ENDFOR
         Ilinha_beta_lambda=T_beta_k##TRANSPOSE(E_lambda_k_gama)
      ENDIF

      ;method 2 for the feature suppression and enhancement process
      IF (feature_method EQ '2') THEN BEGIN
         PRINTF,2,'Method used for the feature suppression and enhancement procedure: '+feature_method
         PRINTF,2,''
         gama_k_list=''
         eigenvectors_list=''
         eigenvalues_list=''
         READ,eigenvectors_number,PROMPT='Enter the number of the eigenvectors to be considered in the feature suppression and enhancement process: '
         READ,eigenvectors_list,PROMPT='Enter the name of the file containing the eigenvectors to be considered in the process: '
         READ,gama_k_list,PROMPT='Enter the name of the file containing the feature factors to be attributed to the eigenvectors in the process: '
         READ,eigenvalues_list,PROMPT='Enter the name of the file containing the eigenvalues of the eigenvectors to be considered in the feature enhancement process: '
         gama_k=make_array(eigenvectors_number,/FLOAT)
         eigenvectors=make_array(eigenvectors_number,/INTEGER)
         eigenvalues=make_array(eigenvectors_number,/FLOAT)
         OPENR,1,infolder+gama_k_list
         READF,1,gama_k
         close,1
         OPENR,1,infolder+eigenvectors_list
         READF,1,eigenvectors
         close,1
         OPENR,1,infolder+eigenvalues_list
         READF,1,eigenvalues
         close,1
         PRINTF,2,'Eigenvectors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,eigenvectors[w]
         ENDFOR
         PRINTF,2,''
         PRINTF,2,'Feature factors associated to each one of the above eigenvectors: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,gama_k[w]
         ENDFOR
         PRINTF,2,''
         E_lambda_k=make_array(eigenvectors_number,m,/FLOAT)
         T_beta_k=make_array(eigenvectors_number,n,/FLOAT)
         E_lambda_k_gama=make_array(eigenvectors_number,m,/FLOAT)
         T_beta_k_N=make_array(eigenvectors_number,n,/FLOAT)
         eigenspectrum=make_array(2,m,/FLOAT)
         tomogram=make_array(maxspatpx_x-minspatpx_x+1,maxspatpx_y-minspatpx_y+1,/FLOAT)
         FOR w=0, eigenvectors_number-1 DO BEGIN
            l=STRING(eigenvectors[w])
            l=STRTRIM(l,1)
            OPENR,1,infolder+eigenspectra_prefix+l+'.txt'
            READF,1,eigenspectrum
            close,1
            E_lambda_k[w,*]=eigenspectrum[1,*]
            tomogram=MRDFITS(infolder+tomograms_prefix+l+'.fits')
            y=0
            FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
               FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
                  T_beta_k[w,y]=tomogram[i,j]
                  y=y+1
               ENDFOR
            ENDFOR
         ENDFOR
         FOR w=0, eigenvectors_number-1 DO BEGIN
            E_lambda_k_gama[w,*]=E_lambda_k[w,*]*gama_k[w]
         ENDFOR
         PRINTF,2,'Nk factors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            Nk=1/(eigenvalues[w]*(n-1))^(0.5)
            PRINTF,2,Nk
            T_beta_k_N[w,*]=T_beta_k[w,*]*Nk
         ENDFOR
         Ilinha_beta_lambda=T_beta_k_N##TRANSPOSE(E_lambda_k_gama)
      ENDIF
   ENDIF

   ;feature suppression and enhancement executed using the tables
   IF (method EQ 'tables') or (method EQ 't') THEN BEGIN
      SCORE=''
      PC=''
      feature_method=''
      READ,SCORE,PROMPT="Enter the file's name of the table containing the parameters of the tomograms obtained with the PCA tomography: "
      READ,PC,PROMPT="Enter the file's name of the table containing the parameters of the eigenvectors obtained with the PCA tomography: "
      E_lambda_k=make_array(m,m,/FLOAT)
      T_beta_k=make_array(m,n,/FLOAT)
      OPENR,1,infolder+SCORE
      READF,1,T_beta_k
      close,1
      OPENR,1,infolder+PC
      READF,1,E_lambda_k
      close,1
      READ,feature_method,PROMPT='Do you wish to perform the feature suppression and enhancement process using method 1 or method 2? '
      WHILE (feature_method NE '1') and (feature_method NE '2') DO BEGIN
         print, 'ERROR: you must answer "1" or "2"!'
         READ,feature_method,PROMPT='Do you wish to perform the feature suppression and enhancement process using method 1 or method 2?'
      ENDWHILE

      ;método 1 for the feature suppression and enhancement process
      IF (feature_method EQ '1') THEN BEGIN
         PRINTF,2,'Method used for the feature suppression and enhancement procedure: '+feature_method
         PRINTF,2,''
         eigenvectors_list=''
         weights_list=''
         READ,eigenvectors_number,PROMPT='Enter the number of the eigenvectors to be considered in the feature suppression and enhancement process: '
         READ,eigenvectors_list,PROMPT='Enter the name of the file containing the eigenvectors to be considered in the process: '
         READ,weights_list,PROMPT='Enter the name of the file containing the feature factors to be attributed to the eigenvectors in the process: '
         weights=make_array(eigenvectors_number,/FLOAT)
         gama_k=make_array(m,/FLOAT)
         eigenvectors=make_array(eigenvectors_number,/INTEGER)
         OPENR,1,infolder+weights_list
         READF,1,weights
         close,1
         OPENR,1,infolder+eigenvectors_list
         READF,1,eigenvectors
         close,1
         PRINTF,2,'Eigenvectors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,eigenvectors[w]
         ENDFOR
         PRINTF,2,''
         PRINTF,2,'Feature factors associated to each one of the above eigenvectors: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,weights[w]
         ENDFOR
         PRINTF,2,''
         E_lambda_k_gama=make_array(m,m,/FLOAT)
         FOR w=0, m-1 DO BEGIN
            gama_k[w]=0.0
         ENDFOR
         FOR w=0, eigenvectors_number-1 DO BEGIN
            gama_k[eigenvectors[w]-1]=weights[w]
         ENDFOR
         FOR w=0, m-1 DO BEGIN
            E_lambda_k_gama[w,*]=E_lambda_k[w,*]*gama_k[w]
         ENDFOR
         Ilinha_beta_lambda=T_beta_k##TRANSPOSE(E_lambda_k_gama)
      ENDIF

      ;method 2 for the feature suppression and enhancement process
      IF (feature_method EQ '2') THEN BEGIN
         PRINTF,2,'Method used for the feature suppression and enhancement procedure: '+feature_method
         PRINTF,2,''
         eigenvectors_list=''
         weights_list=''
         eigenvalues_list=''
         READ,eigenvectors_number,PROMPT='Enter the number of the eigenvectors to be considered in the feature suppression and enhancement process: '
         READ,eigenvectors_list,PROMPT='Enter the name of the file containing the eigenvectors to be considered in the process: '
         READ,weights_list,PROMPT='Enter the name of the file containing the feature factors to be attributed to the eigenvectors in the process: '
         READ,eigenvalues_list,PROMPT='Enter the name of the file containing the eigenvalues of the eigenvectors to be considered in the feature enhancement process: '
         weights=make_array(eigenvectors_number,/FLOAT)
         gama_k=make_array(m,/FLOAT)
         eigenvectors=make_array(eigenvectors_number,/INTEGER)
         eigenvalues=make_array(eigenvectors_number,/FLOAT)
         OPENR,1,infolder+weights_list
         READF,1,weights
         close,1
         OPENR,1,infolder+eigenvectors_list
         READF,1,eigenvectors
         close,1
         OPENR,1,infolder+eigenvalues_list
         READF,1,eigenvalues
         close,1
         PRINTF,2,'Eigenvectors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,eigenvectors[w]
         ENDFOR
         PRINTF,2,''
         PRINTF,2,'Feature factors associated to each one of the above eigenvectors: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            PRINTF,2,weights[w]
         ENDFOR
         PRINTF,2,''
         E_lambda_k_gama=make_array(m,m,/FLOAT)
         T_beta_k_N=make_array(m,n,/FLOAT)
         FOR w=0, m-1 DO BEGIN
            gama_k[w]=0.0
         ENDFOR
         FOR w=0, eigenvectors_number-1 DO BEGIN
            gama_k[eigenvectors[w]-1]=weights[w]
         ENDFOR
         FOR w=0, m-1 DO BEGIN
            E_lambda_k_gama[w,*]=E_lambda_k[w,*]*gama_k[w]
         ENDFOR
         T_beta_k_N=T_beta_k
         PRINTF,2,'Nk factors used in the process: '
         FOR w=0, eigenvectors_number-1 DO BEGIN
            Nk=1/(eigenvalues[w]*(n-1))^(0.5)
            PRINTF,2,Nk
            T_beta_k_N[eigenvectors[w]-1,*]=T_beta_k[eigenvectors[w]-1,*]*Nk
         ENDFOR
         Ilinha_beta_lambda=T_beta_k_N##TRANSPOSE(E_lambda_k_gama)
      ENDIF
   ENDIF
   y=0
   FOR i=0, maxspatpx_x-minspatpx_x DO BEGIN
      FOR j=0, maxspatpx_y-minspatpx_y DO BEGIN
         final_cube[i,j,*]=Ilinha_beta_lambda[*,y]
         y=y+1
      ENDFOR
   ENDFOR
   MWRFITS,final_cube,finalcube,headerout
ENDIF
close,2
print, 'Done!'

END



