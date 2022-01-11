pro corrige_header

pasta = 'C:\Users\incta\Desktop\share\4151_cubes\images_archives\mask\paper\si_almudena\n4151_velslices\'
arquivo = 'ngc4151_k035_pa0.fits'
;lambdazero = 20068.7177 ; Ver valor do CRVAL3 no header.
;deltalambda =  2.13384032249451 ; Ver valor do CD3_3 ou CDELT3 no header.
;pixel_referencia =1 ; Ver valor do CRPIX3 no header.
   lambdazero=19650
   deltalambda=2.5
   pixel_referencia=1
   unidade=STRING('angstroms')
   unidade3=STRING('angstroms')
   unidade2=STRING('wtype=linear label=Wavelength units=angstroms units_display=angstrom')
   ;unidade5=600.0
   ;unidade6=1.099
   ;lambdazero2=19388.7324
   ;deltalambda2=2.4362014
   ;pixel_referencia2=1.0

   ;unidade1=-3.47222222222222e-05
  ; unidade2=0.0
  ; unidade3=0.0
  ; unidade4=0.0
   ;unidade5=3.47222222222222e-05
  ; unidade6=0.0
  ; unidade7=0.0
  ; unidade8=0.0
  ; unidade9=1.95
;__________________________________________________________________________________________

cubo = mrdfits(pasta+'\'+arquivo,0,header)
fxaddpar, header, 'CRVAL3', lambdazero
fxaddpar, header, 'CDELT3', deltalambda
fxaddpar, header, 'CDELT1', deltalambda
fxaddpar, header, 'CRPIX3', pixel_referencia
fxaddpar, header, 'CUNIT1', STRING('angstroms')
fxaddpar, header, 'CUNIT2', STRING('angstroms')
fxaddpar, header, 'CUNIT3', STRING('angstroms')
fxaddpar, header, 'CD3_3', deltalambda
;sxaddpar, header, 'CUNIT1', unidade
;sxaddpar, header, 'CUNIT2', unidade
;sxaddpar, header, 'CUNIT3', unidade

;fxaddpar, header, 'CD1_1', unidade1
;fxaddpar, header, 'CD1_2', unidade2
;fxaddpar, header, 'CD1_3', unidade3
;fxaddpar, header, 'CD2_1', unidade4
;fxaddpar, header, 'CD2_2', unidade5
;fxaddpar, header, 'CD2_3', unidade6
;fxaddpar, header, 'CD3_1', unidade7
;fxaddpar, header, 'CD3_2', unidade8
;fxaddpar, header, 'CD3_3', unidade9

;fxaddpar, header, 'EXPTIME', 600.0
;fxaddpar, header, 'AIRMASS', 1.099


fxaddpar, header, 'WAT1_001', unidade2
;____________________________________________________________problema cubos sinfoni eso428
;fxaddpar, header, 'CRPIX1', pixel_referencia2    ;tirar no caso de cubos 3D
;fxaddpar, header, 'CDELT1', deltalambda2     ;idem
;fxaddpar, header, 'CRVAL1', lambdazero2    ;idem
;fxaddpar, header, 'CUNIT1', unidade3	           ;idem
;fxaddpar, header, 'CD1_1', unidade5        ;idem

header = header[6:*]

mwrfits, cubo, pasta+'\h'+arquivo, header

end