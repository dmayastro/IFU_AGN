pro reamostragem_cubo

;nome do cubo de dados a ser reamostrado
cube='C:\Users\incta\Desktop\share\4151_cubes\images_archives\mask\paper\si_almudena\n4151_velslices\hngc4151_k035_pa0zg.fits'

;nome do cubo de dados reamostrado a ser gerado
saida='C:\Users\incta\Desktop\share\4151_cubes\images_archives\mask\paper\si_almudena\n4151_velslices\hngc4151_k035_pa0zgr.fits'


px_or_x = 0.035
px_or_y = 0.035
px_saida_x = 0.018
px_saida_y = 0.018




;início do programa
cubo=MRDFITS(cube,0,header)
z=size(cubo)
Nx = (z[1]*px_or_x/px_saida_x)
Ny = (z[2]*px_or_y/px_saida_y)
if Nx - fix(Nx) GT 0.5 then Nx = fix(Nx+1) else Nx = fix(Nx)
if Ny - fix(Ny) GT 0.5 then Ny = fix(Ny+1) else Ny = fix(Ny)
print, z[1]*px_or_x/Nx, z[2]*px_or_y/Ny
result1=make_array(Nx,Ny,z[3],/FLOAT)
result2=make_array(Nx,Ny,z[3],/FLOAT)

FOR k=0, z[3]-1 DO BEGIN
   result1[*,*,k]=FREBIN(cubo[*,*,k], Nx, Ny,/TOTAL)
   print, k+1
ENDFOR

FOR k=0, z[3]-1 DO BEGIN
   FOR j=0, Ny-1 DO BEGIN
      result2[*,j,k]=INTERPOL(result1[*,j,k],Nx,/LSQUADRATIC)
   ENDFOR
   print, k+1
ENDFOR

FOR k=0, z[3]-1 DO BEGIN
   FOR i=0, Nx-1 DO BEGIN
      result2[i,*,k]=INTERPOL(result2[i,*,k],Ny,/LSQUADRATIC)
   ENDFOR
   print, k+1
ENDFOR

headerout = header[6:*]

MWRFITS,result2,saida,headerout,/CREATE

print,'Pixel x = ',z[1]*px_or_x/Nx,'Pixel y = ',z[2]*px_or_y/Ny

END