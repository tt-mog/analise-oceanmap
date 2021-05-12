# Análise Oceanográfica - Mapas Espaciais
Repositório para scripts de mapas para avaliações oceanográficas

[![N|Solid](https://www.tetratech.com/cs/ttcom/img/tetratech-50th.png)](https://tetratech.com)

### Baixa estrutura 
Clique em Code e Download ZIP para baixar arquivo zipado 

```sh
$ git clone https://github.com/tt-mog/analise-oceanmap.git
```
### Cria o Ambiente Virtual do python
```sh
$ cd analise-oceanmap/
Verifique qual python está sendo utilizado (deve-se usar o python3)
$ which python
$ pip install virtualenv
$ virtualenv env
$ source env/bin/activate  # para utilizar o python local (ambiente virtual)
$ pip install -r requirements.txt # Instala os pacotes necessários
#Instalando dependencias 
$ conda install -c conda-forge xarray dask netCDF4 bottleneck 
$ conda install -c conda-forge cartopy

```

### Cria o Ambiente Virtual do python
```sh
Executando script
$ ipython 

Para gerar mapas MERCATOR:
$ run plot_map_mercator.py 0 0 # Onde o primeiro argumento é a profundidade a ser analisada e o segundo o tipo de quiver (exemplos abaixo)

Para gerar mapas WOA18:
$ run plot_map_WOA18 0 # Onde o primeiro argumento é a profundidade a ser analisada

```

### Resultado
![alt text](https://github.com/tt-mog/analise-oceanmap/blob/master/resultados/corr_sup.png)
![alt text](https://github.com/tt-mog/analise-oceanmap/blob/master/resultados/corr_1000m.png)
![alt text](https://github.com/tt-mog/analise-oceanmap/blob/master/resultados/temp_sal_merc.png)
![alt text](https://github.com/tt-mog/analise-oceanmap/blob/master/resultados/temp_sal_woa18.png)



