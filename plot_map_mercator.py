# -*- coding: utf-8 -*-
"""
Rotina para plotar mapa espacial anual em superficie do MERCATOR a partir do threeds
obs: Incluido os pontos (instrumentos) dos fundeios

Projeto MARGEM EQUATORIAL (MEQ - PREMIER)
Renan Pimentel - mai/2020
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import cmocean
import cartopy
import cartopy.crs as crs
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature, NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pandas as pd
from datetime import datetime
from datetime import date, timedelta
from pathlib import Path
import pytz
import sys
pd.options.mode.chained_assignment = None
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
#from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
#                 cartopy_xlim, cartopy_ylim,smooth2d,)
from matplotlib.colors import LinearSegmentedColormap


#_____________________________________________________________#

#                          Definitions                        #

#_____________________________________________________________#

# Definindo uma função para geração dos mapas
def mapa_base():
    ax = plt.axes(projection=crs.PlateCarree())

    # -- Put a background image on for nice sea rendering --
    gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                        linewidth=0.5, color='black', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    gl.xlabel_style = {'size': 8, 'color': 'black'}
	
    gl.xlocator = mticker.FixedLocator(np.arange(lonmin-2,lonmax+2,2).round())
    gl.ylocator = mticker.FixedLocator(np.arange(latmin-2,latmax+2,2).round())


    states = NaturalEarthFeature(category="cultural", scale="50m",
                            facecolor="beige",
                            name="admin_1_states_provinces_shp")

    ax.add_feature(states, linewidth=1.0, edgecolor="black")
    ax.coastlines('50m', linewidth=1.0)

    mapa_amsul = ShapelyFeature(Reader('shapes/brasil_america.shp').geometries(),
                                    crs.PlateCarree(),
                                    edgecolor='black',
                                    facecolor='beige',
                                    linewidth=0.5)
    ax.add_feature(mapa_amsul)

    return ax


def choose_variable(variavel):
    if variavel == 'Corrente':
        output = ['uo','vo']
        levels = np.arange(0,1.8,0.2)
        ticks = np.arange(0,1.8,0.2)
        cmap = plt.cm.Blues
        modelo = 'MERCATOR'
        label = 'Velocidade (m/s)'

    elif variavel == 'Temperatura':
        output = ['thetao']
        levels = np.arange(20,30.2,0.5)
        ticks = np.arange(20,30.2,1)
        cmap = plt.cm.jet
        modelo = 'MERCATOR'
        label = 'Temperatura (°C)'

    elif variavel == 'Salinidade':
        output = ['so']
        levels = np.arange(20,39.5,1.05)
        ticks = np.arange(20,39.1,1)
        cmap = cmocean.cm.deep
        modelo = 'MERCATOR'
        label = 'Salinidade (psu)'


    return output, levels, cmap, modelo, label, ticks

#_____________________________________________________________#


if __name__ == '__main__':


    lonmin = -54.01
    lonmax = -33.99
    latmin = -6.01
    latmax = 10.01


    variaveis = ['Corrente','Temperatura','Salinidade']
    for variavel in variaveis:
        print('__________________')
        print('Processando '+ variavel)
        output, levels, cmap, modelo, label, ticks = choose_variable(variavel)

        # Setando o index da profundidade
        profundidade = int(sys.argv[1])
        
        # Seleciona tipo de quiver
        tipo = int(sys.argv[2])


        if modelo == 'MERCATOR':
            ds = xr.open_dataset('http://10.50.10.64:8080/thredds/dodsC/3D/mercator_3D_annualmean_2016_2020.nc')

        # Seleciona variavel
        ds_temp = ds[output]

        #Seleciona profundidade
        ds_temp = ds_temp.isel(depth = profundidade)
        ds_sel = ds_temp.squeeze()
        profundidade = int(ds_sel.depth.values)

        # Fazendo calculos necessarios
        if variavel == 'Corrente':
            if (profundidade == 0):
                tamanho = 15 #tamanho de escala
                vscale = 0.5 # quiver de referencia
                escala = 8 # reamostragem

            else:
                tamanho = 10
                vscale = 0.2 
                escala = 4
                levels = np.arange(0,0.6,0.1)
                ticks = np.arange(0,0.6,0.1)
            
            titulo = 'Corrente em {}m \n'.format(str(profundidade),modelo) 

            ds_estacao = xr.Dataset()
            ds_sel['vel'] = (ds_sel.uo**2+ds_sel.vo**2)**0.5
            ds_estacao['vel'] = ds_sel.vel
            ds_estacao['UC'] = ds_sel.uo
            ds_estacao['VC'] = ds_sel.vo
      
        elif variavel == 'Temperatura':
            titulo = 'Temperatura em {}m \n{}'.format(str(profundidade),modelo)  

            ds_estacao = xr.Dataset()
            ds_estacao['temp'] = ds_sel.thetao
     
        elif variavel == 'Salinidade':
            titulo = 'Salinidade em {}m \n{}'.format(str(profundidade),modelo)  

            ds_estacao = xr.Dataset()
            ds_estacao['sal'] = ds_sel.so
        
        X, Y = np.meshgrid(ds_estacao.longitude, ds_estacao.latitude)
        ds_estacao=np.squeeze(ds_estacao)
    
        # Importando base de shapefile
        ax = mapa_base()


        # --------- PLOTANDO BATIMETRIA ---------- #

        import shapefile as shp
        niveis_plotados=[200,1000,2000,3000]
        sf = shp.Reader('shapes/linhas_de_batimetria_meq.shp')

        for shape in sf.shapeRecords():
            if float(shape.record['depth']) in niveis_plotados:
                x = [i[0] for i in shape.shape.points[:]]
                y = [i[1] for i in shape.shape.points[:]]
                ax.plot(x,y,color='black', linestyle='solid',linewidth=0.6)

        # ----- Importando fundeios ---- #
        # ---- POT ---- #
        ax.plot(-37.049,-3.975,'o',markersize=5,transform=crs.PlateCarree(),label = 'P7',color='yellow')
        ax.plot(-37.209,-4.114,'o',markersize=5,transform=crs.PlateCarree(),label = 'P6',color='yellow')

        # ---- CE ---- #
        ax.plot(-39.0078,-2.7527,'o',markersize=5,transform=crs.PlateCarree(),label = 'NOVO CE2',color='blue')
        ax.plot(-38.9319,-2.6190,'o',markersize=5,transform=crs.PlateCarree(),label = 'NOVO CE3',color='blue')
        ax.plot(-39.0478,-2.8056,'o',markersize=5,transform=crs.PlateCarree(),label = 'NOVO CE1',color='blue')
        ax.plot(-39.197,-2.603,'o',markersize=5,transform=crs.PlateCarree(),label = 'WS CE',color='blue')

        if variavel == 'Corrente':

            id_lat = np.where(Y[:,1]>-4)[0][0]
            id_lon = np.where(X[1,:]>-52)[0][0]
            qleg = (id_lat,id_lon)
            qleg_tx = (id_lat+2,id_lon)

            temp = plt.contourf(X,Y,ds_estacao['vel'],cmap=cmap,
                                levels=levels, transform=crs.PlateCarree(), extend='max')
            norm_u = ds_estacao['UC']/np.sqrt(ds_estacao['UC']**2 + ds_estacao['VC']**2)
            norm_v = ds_estacao['VC']/np.sqrt(ds_estacao['UC']**2 + ds_estacao['VC']**2)

            if (tipo == 0):
                # Normalize the arrows:
                quiv = plt.quiver(X[::escala, ::escala],Y[::escala, ::escala],
                norm_u[::escala, ::escala],norm_v[::escala, ::escala],color='k',headwidth=2,headlength=3,units='width',width=0.0015,scale=tamanho)
                clim = (-.02,0.2)
                quiv.set_clim(clim)
                plt.quiver(X[qleg],Y[qleg],vscale,0,color='k',units='width',width=0.0015,scale=tamanho,zorder=50000)
                plt.text(X[qleg_tx],Y[qleg_tx],'%s m/s'%vscale,fontsize=8)
            else:
                tamanho = 2.5
                escala = 5
                quiv = plt.quiver(X[::escala, ::escala],Y[::escala, ::escala],
                norm_u[::escala, ::escala],norm_v[::escala, ::escala],units='xy', angles='xy', scale_units='xy', scale = tamanho, color='black')

        elif variavel == 'Temperatura':
            temp = plt.contourf(X,Y,ds_estacao['temp'],cmap=cmap,
                                    levels=levels, transform=crs.PlateCarree(), extend='both')
 
        elif variavel == 'Salinidade':
            temp = plt.contourf(X,Y,ds_estacao['sal'],cmap=cmap,
                                    levels=levels, transform=crs.PlateCarree(), extend='both')
            

        cbar=plt.colorbar(temp, ax=ax, orientation="vertical", pad=.05)
        cbar.set_label(label,fontsize=10)
        cbar.set_ticks(ticks)

        ax.set_extent([lonmin, lonmax, latmin, latmax])

        plt.title('{} \n Média Anual'.format(titulo), fontsize=10)
    
        from matplotlib.image import imread
        img = imread('Tt_corp_logo_horz_blu.png')
        ax.imshow(img, origin='upper', transform=crs.PlateCarree(), extent=[lonmin,lonmin*0.915, latmin*1.06, latmin*0.7],zorder=10)

        if modelo == 'MERCATOR':
            if variavel == 'Corrente':
                if (tipo == 0):
                    plt.savefig('resultados/mapa_mercator_{}_{}m_annual_quiv0.png'.format(variavel.lower(),str(profundidade)), dpi=500, bbox_inches='tight', pad_inches=0.1)
                else:
                    plt.savefig('resultados/mapa_mercator_{}_{}m_annual_quiv1.png'.format(variavel.lower(),str(profundidade)), dpi=500, bbox_inches='tight', pad_inches=0.1)
            else:
                plt.savefig('resultados/mapa_mercator_{}_{}m_annual.png'.format(variavel.lower(),str(profundidade)), dpi=500, bbox_inches='tight', pad_inches=0.1)

        plt.close()


print('Fim do Script!')

