import pandas as pd
import numpy as np
import sklearn as sk
import matplotlib.pyplot as plt

syntenydata = pd.read_csv('data/pycharm_deneme_dataset.csv')

cldata = syntenydata.drop(['annotation','num','orf_name'], axis=1).fillna(0)

from scipy.cluster import hierarchy
Z=hierarchy.linkage(cldata,method='ward')

plt.figure()
dn = hierarchy.dendrogram(Z)
plt.show()
# from sklearn.cluster import AgglomerativeClustering
#
# clustering = AgglomerativeClustering().fit(cldata)
# clustering




#clustering.labels_