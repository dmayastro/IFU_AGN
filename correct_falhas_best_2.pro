pro correct_falhas_best_2

;cubo a ser corrigido
entrada='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zc.fits'

;cubo corrigido a ser gerado
saida='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zcg.fits'

;arquivo contendo os parâmetros de todas as correções. Esse arquivo deve ser configurado de tal forma que
;cada linha contenha os parâmetros correspondentes a cada uma das correções a serem feitas, ou seja, o número de linhas
;do arquivo será igual ao número de correções. Cada uma dessas linhas do arquivo deve conter os parâmetros na seguinte ordem:
;Ni (pixel espectral inicial do intervalo no qual será feita a interpolação), Nf (pixel espectral final do intervalo no qual
;será feita a interpolação), Nxi (pixel espacial inicial no eixo x a ser considerado na interpolação), Nxf (pixel espacial final
;no eixo x a ser considerado na interpolação), Nyi (pixel espacial inicial no eixo y a ser considerado na interpolação), Nyf (pixel
;espacial final no eixo y a ser considerado na interpolação), Nmediai (número de pixeis espectrais antes de Ni a serem considerados
;no cálculo do valor médio que será utilizado para realizar a interpolação), Nmediaf (número de espectrais após Nf a serem
;considerados no cálculo do valor médio que será utilizado para realizar a interpolação)
lista_correcoes='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\gap.txt'

;número de correções a serem feitas (que é igual ao número de linhas do arquivo dado pelo parâmetro anterior)
numero_correcoes=7


;início do programa

image=MRDFITS(entrada,0,header)
correcoes=make_array(8,numero_correcoes,/FLOAT)
OPENR,1,lista_correcoes
READF,1,correcoes
close,1
z=size(image)
FOR w=0, numero_correcoes-1 DO BEGIN
;loop para cada correção
   Ni=correcoes[0,w]
   Nf=correcoes[1,w]
   Nxi=correcoes[2,w]
   Nxf=correcoes[3,w]
   Nyi=correcoes[4,w]
   Nyf=correcoes[5,w]
   ;aqui Nmedia são os valores de X e Y em uma imagem, para interpolar.
   Nmediai=correcoes[6,w]
   Nmediaf=correcoes[7,w]
   falhain=make_array(Nmediai,/float)
   falhafin=make_array(Nmediaf,/float)
   falha=make_array(2,/float)
   falhacor=make_array(Nf-Ni+1,/float)
   FOR i=Nxi-1, Nxf-1 DO BEGIN
      FOR j=Nyi-1, Nyf-1 DO BEGIN
         FOR l=0, Nmediai-1 DO BEGIN
            falhain[l]=image[i,j,Ni-1-l]
           ;falhain[l]=image[i,j,Ni-1-l]
         ENDFOR
         FOR l=0, Nmediaf-1 DO BEGIN
            falhafin[l]=image[i,j,Nf-1+l]
         ENDFOR
         falha[0]=MEAN(falhain)
         falha[1]=MEAN(falhafin)
         falhacor=INTERPOL(falha,Nf-Ni+1)
         FOR l=0, Nf-Ni DO BEGIN
           image[i,j,Ni-1+l]=falhacor[l]
           ; image[i,j,Ni-1+l]=0.0
         ENDFOR
      ENDFOR
   ENDFOR
   print, w
ENDFOR

MWRFITS,image,saida,header,/CREATE

END
