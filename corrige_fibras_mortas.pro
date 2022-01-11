pro corrige_fibras_mortas

cube='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\std\std_ESO\nbcmsifs20170808.098.cube.fits'
mask='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\std\std_ESO\nbcmsifs20170808.096.cube_mask.fits'
saida='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\std\std_ESO\nbcmsifs20170808.098.cubec.fits'

cubo=MRDFITS(cube,0,header_cubo)
z=size(cubo)
mascara=MRDFITS(mask,0,header)
cubo_corrigido=cubo
limites_falha=make_array(2,/FLOAT)

FOR k=0, z[3]-1 DO BEGIN
   mascara_corrigida=mascara
   FOR i=2, z[1]-2 DO BEGIN
      FOR j=0, z[2]-1 DO BEGIN
         IF (mascara_corrigida[i,j] EQ 1.0) THEN BEGIN
            x_menos=i-1
            WHILE (mascara_corrigida[x_menos,j] EQ 1.0) DO BEGIN
               x_menos=x_menos-1
            ENDWHILE
            x_mais=i+1
            WHILE (mascara_corrigida[x_mais,j] EQ 1.0) DO BEGIN
               x_mais=x_mais+1
            ENDWHILE
            limites_falha[0]=cubo_corrigido[x_menos,j,k]
            limites_falha[1]=cubo_corrigido[x_mais,j,k]
            cubo_corrigido[x_menos:x_mais,j,k]=INTERPOL(limites_falha,x_mais-x_menos+1)
            mascara_corrigida[x_menos:x_mais,j]=0.0
         ENDIF
      ENDFOR
   ENDFOR
ENDFOR

MWRFITS,cubo_corrigido,saida,header_cubo,/CREATE

END






