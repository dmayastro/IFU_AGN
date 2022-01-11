function determine_goodPixels, logLam, lamRangeTemp, redshift

lines = [21218d,22232d,21661d,22478d] ; Wavelengths of Hbeta and [OIII] doublet
dv = [300d,300d,250d,300d] ; width/2 of masked gas emission region in km/s
c = 299792.458d ; speed of light in km/s

flag = bytarr(n_elements(logLam))

for j=0,n_elements(lines)-1 do $
    flag or= logLam gt alog(lines[j]) + (redshift - dv[j])/c $
         and logLam lt alog(lines[j]) + (redshift + dv[j])/c

flag or= logLam lt alog(lamRangeTemp[0]) + (redshift + 900d)/c ; Mask edges of
flag or= logLam gt alog(lamRangeTemp[1]) + (redshift - 900d)/c ; stellar library

return, where(flag eq 0)
end

pro obtem_mapas_incertezas_simulacao_new

cubo='/sto/home/robertobm/NGC6951_incertezas_teste2/MrhE19rzhc_ppxf80.fits'
Ni_SdivN=820
Nf_SdivN=870
diretorio_templates='/sto/home/robertobm/NIFS_rebin_original/'
objeto='NGC6951'
outdir='/sto/home/robertobm/NGC6951_incertezas_teste2/'

cube=MRDFITS(cubo,0,h)
z=size(cube)
galaxy_redshift=0d
initial_guess_dispersion=50d
spectrum=make_array(z[3],/DOUBLE)
velocidades=make_array(20,/FLOAT)
dispersoes=make_array(20,/FLOAT)
h3s=make_array(20,/FLOAT)
h4s=make_array(20,/FLOAT)
lambdas_SdivN=make_array(Nf_SdivN-Ni_SdivN+1.0,/FLOAT)
lambdazero=FXPAR(h,'CRVAL3')
deltalambda=FXPAR(h,'CDELT3')
reference_pixel=FXPAR(h,'CRPIX3')

FOR n=Ni_SdivN, Nf_SdivN DO BEGIN
   lambdas_SdivN[n-Ni_SdivN]=lambdazero+(n-reference_pixel)*deltalambda
ENDFOR

lamRange1 = sxpar(h,'CRVAL3') + [0d,sxpar(h,'CDELT3')*(sxpar(h,'NAXIS3')-1d)]
spectrum[*]=cube[0,0,*]
log_rebin, lamRange1, spectrum, galaxy, logLam1, VELSCALE=velScale
galaxy = galaxy/median(galaxy) ; Normalize spectrum to avoid numerical issues
noise = galaxy*0 + 1           ; Assume constant noise per pixel here

vazdekis = file_search(diretorio_templates+"*.fits",COUNT=nfiles)

fits_read, vazdekis[0], ssp, h
lamRange2 = sxpar(h,'CRVAL1') + [0d,sxpar(h,'CDELT1')*(sxpar(h,'NAXIS1')-1d)]
log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
templates = dblarr(n_elements(sspNew),nfiles)

for j=0,nfiles-1 do begin
   fits_read, vazdekis[j], ssp
   ;ssp = convol(ssp,lsf) ; Degrade template to SAURON resolution
   log_rebin, lamRange2, ssp, sspNew, VELSCALE=velScale
   ;templates[*,j] = sspNew/median(sspNew) ; Normalizes templates
   templates[*,j]=sspNew
endfor

c = 299792.458d
dv = (logLam2[0]-logLam1[0])*c

redshift = galaxy_redshift ; Initial estimate of the galaxy redshift in km/s
goodPixels = determine_goodPixels(logLam1,lamRange2,redshift)
w=1

FOR i=0, z[1]-1 DO BEGIN
   FOR j=0, z[2]-1 DO BEGIN
      ;reta=LINFIT(lambdas_SdivN,cube[i,j,Ni_SdivN-1:Nf_SdivN-1])
      ;linear_spectrum=reta[0]+reta[1]*lambdas_SdivN
      ;PLOT, lambdas_SdivN, linear_spectrum
      ;OPLOT, lambdas_SdivN, cube[i,j,Ni_SdivN-1:Nf_SdivN-1]

      ;sum=0.0
      ;FOR k=Ni_SdivN-1, Nf_SdivN-1 DO BEGIN
      ;   sum=sum+(cube[i,j,k]-linear_spectrum[k+1-Ni_SdivN])^2.0
      ;ENDFOR
      ;rms=SQRT(sum/(Nf_SdivN-Ni_SdivN+1.0))
      ;print, rms
      spectrum[*]=cube[i,j,*]
      log_rebin, lamRange1, spectrum, galaxy, logLam1, VELSCALE=velScale
      galaxy = galaxy/median(galaxy) ; Normalize spectrum to avoid numerical issues
      noise = galaxy*0 + 1           ; Assum
      start = [redshift, initial_guess_dispersion] ; (km/s), starting guess for [V,sigma]
      ppxf, templates, galaxy, noise, velScale, start, sol, BESTFIT=bestFit, $
      GOODPIXELS=goodPixels, MOMENTS=4, DEGREE=4, $
      VSYST=dv, ERROR=error
      velocidade=sol[0]
      dispersao=sol[1]
      h3=sol[2]
      h4=sol[3]
      residuals = galaxy-bestFit
      histograma=HISTOGRAM(residuals,NBINS=50,LOCATIONS=locations)
      Result = GAUSSFIT(locations, histograma, A)
      ;print, A[2]
      ;seed=noise
      FOR n=0, 19 DO BEGIN
         new_noise=DOUBLE(RANDOMN(n+1,z[3])*A[2])
         galaxy=bestFit+new_noise
         ppxf, templates, galaxy, noise, velScale, start, sol, $
         GOODPIXELS=goodPixels, MOMENTS=4, DEGREE=4, $
         VSYST=dv, ERROR=error
         velocidades[n]=sol[0]
         dispersoes[n]=sol[1]
         h3s[n]=sol[2]
         h4s[n]=sol[3]
         seed=noise
      ENDFOR
      ;residuos_velocidades=velocidades-velocidade
      ;residuos_dispersoes=dispersoes-dispersao
      ;residuos_h3s=h3s-h3
      ;residuos_h4s=h4s-h4
      s=string(w)
      s=strtrim(s,1)
      OPENW,1,outdir+'mapa_incertezas_'+objeto+'_'+s+'.txt'
      PRINTF,1,STDDEV(velocidades)
      PRINTF,1,STDDEV(dispersoes)
      PRINTF,1,STDDEV(h3s)
      PRINTF,1,STDDEV(h4s)
      close,1
      w=w+1
   ENDFOR
ENDFOR


END
