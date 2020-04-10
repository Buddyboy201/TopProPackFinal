import matplotlib as plt
import seaborn as sns
from collections import Counter
import centroid_protein
import streamlit
import matplotlib.pyplot as plt
import pandas as pd


def draw_countplot(data, x, title):
	df = pd.DataFrame()
	df[x] = data
	ax = sns.countplot(x=x, data=df, order=list(range(max(data)+1)))
	plt.title("Clique Sizes Histogram")
	for p in ax.patches:
		ax.annotate('{}'.format(p.get_height()), (p.get_x() + 0.1, p.get_height()))
	plt.show()

def draw_histogram(data, title, normalized=False):
	ax = None
	if not normalized:
		ax = sns.distplot(data, kde=False, norm_hist=False)
	else:
		ax = sns.distplot(data)
	plt.show()

def draw_heatmap(heatmap_data):
	sns.heatmap(heatmap_data)
	plt.show()

#draw_countplot([len(s) for s in centroid_protein.P.centroid_cliques], "clique_sizes", "Clique Size Frequencies")
#draw_histogram(centroid_protein.P.centroid_clique_distances, "Centroid Clique Distances", normalized=True)
