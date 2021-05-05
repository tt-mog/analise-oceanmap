# -*- coding: utf-8 -*-
"""
Rotina para ler o Netcdf do WOA e plotar mapa com pontos
dos fundeios e batimetria

Projeto MARGEM EQUATORIAL (MEQ - PREMIER)
Renan Pimentel - maio/2020
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
import cartopy.feature as cfeature
import pandas as pd
from datetime import datetime
from datetime import date, timedelta
from pathlib import Path
#from PIL import Image
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

    gl.xlabels_top = False
    gl.ylabels_right = False
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
    #ax.add_feature(mapa_amsul)

    mapa_amsul = ShapelyFeature(Reader('shapes/brasil_america.shp').geometries(),
                                    crs.PlateCarree(),
                                    edgecolor='black',
                                    facecolor='beige',
                                    linewidth=0.5)
    ax.add_feature(mapa_amsul)

    return ax

#_____________________________________________________________#


if __name__ == '__main__':

    lonmin = -54.01
    lonmax = -33.99
    latmin = -6.01
    latmax = 10.01


    if variavel == 'TSM':
        ds = xr.open_dataset('../DADOS/WOA_MEq_TEMP_completo.nc',decode_times=False)

        output = ['TEMP']
        #cmap = cmocean.cm.thermal
        cmap = plt.cm.RdBu_r
        modelo = 'WOA18'
        label = 'Temperature (°C)'

        # Seleciona variavel
        ds_sel = ds[output]
        ds_sel = ds_sel.isel(DEPTH=0)
        ds_sel = ds_sel.squeeze()
        
        ds_estacao=np.squeeze(ds_sel)
        X, Y = np.meshgrid(ds_estacao.LON, ds_estacao.LAT)

        # Importando base de shapefile
        ax = mapa_base()

        levels = np.arange(20,30.2,0.2)
        
        temp = plt.contourf(X,Y,ds_estacao['TEMP'],cmap=cmap,
        levels=levels, transform=crs.PlateCarree(), extend='both')
        
        cbar=plt.colorbar(temp, ax=ax, orientation="vertical", pad=.05)
        cbar.set_label(label,fontsize=10)
        cbar.set_ticks(np.arange(20,32,2))
        ax.set_extent([lonmin, lonmax, latmin, latmax])
        plt.title('Temperatura Climatológica - NODC / WOA18')


    elif variavel == 'Salinidade':
        ds = xr.open_dataset('../DADOS/WOA_MEq_SAL_completo.nc',decode_times=False)

        output = ['SAL']
        cmap = cmocean.cm.deep
        modelo = 'WOA18'
        label = 'Salinity (PSU)'

        # Seleciona variavel
        ds_sel = ds[output]
        ds_sel = ds_sel.isel(DEPTH=0)
        ds_sel = ds_sel.squeeze()
        
        ds_estacao=np.squeeze(ds_sel)
        X, Y = np.meshgrid(ds_estacao.LON, ds_estacao.LAT)

        # Importando base de shapefile
        ax = mapa_base()

        levels = np.arange(20,39.5,1.05)
        
        temp = plt.contourf(X,Y,ds_estacao['SAL'],cmap=cmap,
        levels=levels, transform=crs.PlateCarree(), extend='both')
        
        cbar=plt.colorbar(temp, ax=ax, orientation="vertical", pad=.05)
        cbar.set_label(label,fontsize=10)
        cbar.set_ticks(np.arange(20,39.1,1))
        ax.set_extent([lonmin, lonmax, latmin, latmax])
        plt.title('Salinidade Climatológica - NODC / WOA18')

    # --------- PLOTANDO BATIMETRIA ---------- #

    import shapefile as shp
    niveis_plotados=[200,1000,2000,3000]
    sf = shp.Reader('../../01_shapes/linhas_de_batimetria_meq.shp')
    #plt.figure()
    for shape in sf.shapeRecords():
        if float(shape.record['depth']) in niveis_plotados:
            x = [i[0] for i in shape.shape.points[:]]
            y = [i[1] for i in shape.shape.points[:]]
            ax.plot(x,y,color='black', linestyle='solid',linewidth=0.6)

    from matplotlib.image import imread
    img = imread('Tt_corp_logo_horz_blu.png')
    ax.imshow(img, origin='upper', transform=crs.PlateCarree(), extent=[lonmin,lonmin*0.915, latmin*1.06, latmin*0.7],zorder=10)
    
    if variavel == 'TSM':
        plt.savefig('../RESULTADOS/MAPAS_Climatologicos/ANUAL/TSM_Climatologico.png', dpi=500, bbox_inches='tight')    
    elif variavel == 'Salinidade':    
        plt.savefig('../RESULTADOS/MAPAS_Climatologicos/ANUAL/Salinity_Climatologico.png', dpi=500, bbox_inches='tight')

    plt.close()
