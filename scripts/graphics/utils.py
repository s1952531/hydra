#!/usr/bin/env python3

import os
homedir=os.getenv('HOME')
rootdir=os.path.join(homedir,'hydra','scripts')
graphicsdir=os.path.join(rootdir,'graphics')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs

def get_colourmap():
    cmap_file = open(os.path.join(graphicsdir,'colourmap'),'r')
    cmap_list = []
    n=254
    for i in range(n):
        file_line = cmap_file.readline().split()
        a,b,c = file_line[0],file_line[1],file_line[2]
        ele = (float(a),float(b),float(c))
        cmap_list.append(ele)
    stamap = clrs.LinearSegmentedColormap.from_list('amap',cmap_list,N=n)
    #plt.register_cmap(cmap=stamap)
    mpl.colormaps.get_cmap(stamap)

    cmap_file = open(os.path.join(graphicsdir,'wbgyr'),'r')
    cmap_list = []    
    for i in range(n):
        file_line = cmap_file.readline().split()
        a,b,c = file_line[0],file_line[1],file_line[2]
        ele = (float(a),float(b),float(c))
        cmap_list.append(ele)
    wbgyr = clrs.LinearSegmentedColormap.from_list('wbgyr',cmap_list,N=n)
    #plt.register_cmap(cmap=wbgyr)
    mpl.colormaps.get_cmap(wbgyr)
    return {'stamap':stamap,'wbgyr':wbgyr,'jet':plt.cm.jet,'seismic':plt.cm.seismic,'rjet':plt.cm.jet_r,
            'rseismic':plt.cm.seismic_r,'autumn':plt.cm.autumn,'winter':plt.cm.winter,'spring':plt.cm.spring,
            'summer':plt.cm.summer,'hot':plt.cm.hot,'bone':plt.cm.bone,'rbone':plt.cm.bone_r,'cool':plt.cm.cool,
            'copper':plt.cm.copper,'gray':plt.cm.gray,'hsv':plt.cm.hsv,'bwr':plt.cm.bwr,'prism':plt.cm.prism,
            'pink':plt.cm.pink,'RdBu':plt.cm.RdBu,'rRdBu':plt.cm.RdBu_r,'Blues':plt.cm.Blues,'rBlues':plt.cm.Blues_r,
            'GnBu':plt.cm.GnBu,'PuBu':plt.cm.PuBu,'BuPu':plt.cm.BuPu,'YlGnBu':plt.cm.YlGnBu,'Greys':plt.cm.Greys,
            'binary':plt.cm.binary,'rainbow':plt.cm.gist_rainbow,'flag':plt.cm.flag}
