pro convoluir_imagem


pasta = 'C:\Users\incta\Desktop\share\eso428g14\calib\extra_maps\
mapa_in = 'HST_WFPC2_Halpha_i
psfmap = 'psf_sum_halpha
saida = 'HST_WFPC2_Halpha_concolved'


;____________________________________________________________________________________________

in_map = mrdfits(pasta+mapa_in+'.fits',0,header)
psf = mrdfits(pasta+psfmap+'.fits',0,header)

tam_ini = size(in_map)

unconvolved_map = in_map

tam = size(unconvolved_map)

convolved_map = convolve(unconvolved_map, psf)

out_map = convolved_map

mwrfits, out_map, pasta+saida+'.fits', /create
end