pro correct_falhas_best_2

;cubo a ser corrigido
entrada='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zc.fits'

;cubo corrigido a ser gerado
saida='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\calib_ngc4507_new_e15zcg.fits'

;arquivo contendo os par�metros de todas as corre��es. Esse arquivo deve ser configurado de tal forma que
;cada linha contenha os par�metros correspondentes a cada uma das corre��es a serem feitas, ou seja, o n�mero de linhas
;do arquivo ser� igual ao n�mero de corre��es. Cada uma dessas linhas do arquivo deve conter os par�metros na seguinte ordem:
;Ni (pixel espectral inicial do intervalo no qual ser� feita a interpola��o), Nf (pixel espectral final do intervalo no qual
;ser� feita a interpola��o), Nxi (pixel espacial inicial no eixo x a ser considerado na interpola��o), Nxf (pixel espacial final
;no eixo x a ser considerado na interpola��o), Nyi (pixel espacial inicial no eixo y a ser considerado na interpola��o), Nyf (pixel
;espacial final no eixo y a ser considerado na interpola��o), Nmediai (n�mero de pixeis espectrais antes de Ni a serem considerados
;no c�lculo do valor m�dio que ser� utilizado para realizar a interpola��o), Nmediaf (n�mero de espectrais ap�s Nf a serem
;considerados no c�lculo do valor m�dio que ser� utilizado para realizar a interpola��o)
lista_correcoes='C:\Users\incta\Desktop\new_sinfoni_nifs_sifs\Alberto\gap.txt'

;n�mero de corre��es a serem feitas (que � igual ao n�mero de linhas do arquivo dado pelo par�metro anterior)
numero_correcoes=7


;in�cio do programa

image=MRDFITS(entrada,0,header)
correcoes=make_array(8,numero_correcoes,/FLOAT)
OPENR,1,lista_correcoes
READF,1,correcoes
close,1
z=size(image)
FOR w=0, numero_correcoes-1 DO BEGIN
;loop para cada corre��o
   Ni=correcoes[0,w]
   Nf=correcoes[1,w]
   Nxi=correcoes[2,w]
   Nxf=correcoes[3,w]
   Nyi=correcoes[4,w]
   Nyf=correcoes[5,w]
   ;aqui Nmedia s�o os valores de X e Y em uma imagem, para interpolar.
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
