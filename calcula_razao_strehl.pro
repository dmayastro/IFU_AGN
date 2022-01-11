pro calcula_razao_Strehl

imagem='C:\Users\dmay\Desktop\shared\ngc7582_steiner\cubes_250\pca_psf_lucy\psf.fits'
Xmin=38.0
Xmax=58.0
Ymin=38.0
Ymax=58.0
;cube='cubofinal_GQLUPUS_reamostrado_cor_filtered1_deconvolved_report_122.fits'

image=MRDFITS(imagem)
image_cut=make_array(Xmax-Xmin+1,Ymax-Ymin+1)
image_cut_norm=make_array(Xmax-Xmin+1,Ymax-Ymin+1)
image_Airy=make_array(Xmax-Xmin+1,Ymax-Ymin+1)
image_Airy_norm=make_array(Xmax-Xmin+1,Ymax-Ymin+1)
image_cut=image[Xmin-1:Xmax-1,Ymin-1:Ymax-1]
image_cut_norm=image_cut/TOTAL(image_cut)
image_Airy=psf_airy_v3(Xmax-Xmin+1, Ymax-Ymin+1, 21661.0e-10, 0.0625, 8.2)
image_Airy_norm=image_Airy/TOTAL(image_Airy)

razao_Strehl=MAX(image_cut_norm)/MAX(image_Airy_norm)

print, razao_Strehl

END

