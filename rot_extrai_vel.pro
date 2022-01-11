pro rot_extrai_vel

imagem='C:\Users\incta\Desktop\share\4151_cubes\4151_banda-k\h2_vel_M6zz_offsetzz.fits'
;x é 0 e y é -90 - rotaciona do topo (N do centro defindo abaixo) para o sentido horário e a medida do pseudo-slit começa a medir da direita para a esquerda.
;PA=106.0 ;eso428
PA=-43
px=0.021
y1=58.0
xc=60.0
yc=58.0
saida='C:\Users\incta\Desktop\share\4151_cubes\4151_banda-k\vel_h2.txt'
saida_imagem='C:\Users\incta\Desktop\share\4151_cubes\4151_banda-k\vel_h2.fits'

angle=PA
image=MRDFITS(imagem,0,header)
image_rotated=ROT(image, angle, 1.0, xc-1, yc-1, INTERP=1, CUBIC=-0.5)
MWRFITS,image_rotated,saida_imagem,/CREATE


z=size(image_rotated)
perfil1=make_array(z[1],/FLOAT)
perfil_medio=make_array(2,z[1],/FLOAT)

perfil_medio[0,*]=(1.0+FINDGEN(z[1]))*px-(z[1]*px/2)

perfil1=image_rotated[*,y1-1]
perfil_medio[1,*]=perfil1

OPENW,1,saida
PRINTF,1,perfil_medio
close,1

END