#load package
import scvelo as scv

#load data
adata = scv.read("mouse2.h5ad")

#load positions, clusters, coloros from SeuratMab21l2
X_tsne = scv.load('cell_embeddings.csv', index_col=0)
adata.obsm['X_tsne'] = X_tsne.values
cluster = scv.load('clusters.csv', index_col=0)
adata.obs['seurat_clusters'] = cluster.values
seurat_colors = scv.load('seurat_colors.csv', index_col=0)
new_colors = ['#F8766D','#B79F00', '#00BA38', '#00BFC4','#619CFF','#F564E3']
adata.uns['seurat_clusters_colors']=new_colors

#AnnData object with n_obs × n_vars = 462 × 54752
#   obs: 'orig.ident', 'nCount_spliced', 'nFeature_spliced', 'nCount_unspliced', 'nFeature_unspliced', 'nCount_ambiguous', 
#        'nFeature_ambiguous', 'nCount_spanning', 'nFeature_spanning', 'nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.0.8', 'seurat_clusters'
#   var: 'features', 'ambiguous_features', 'spanning_features', 'spliced_features', 'unspliced_features'
#   obsm: 'X_tsne'
#   layers: 'ambiguous', 'spanning', 'spliced', 'unspliced'

#preprocess data
#Number of highly-variable genes to keep = 2000
scv.pp.filter_and_normalize(adata, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=20, n_neighbors=20)

#estimate rna velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

#project the velocities
#the velocities are projected onto any embedding, specified by basis, and visualized in one of these ways:
#on cellular level with scv.pl.velocity_embedding, as gridlines with scv.pl.velocity_embedding_grid,
#or as streamlines with scv.pl.velocity_embedding_stream
scv.pl.velocity_embedding(adata, basis="tsne", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
scv.pl.velocity_embedding_grid(adata, basis="tsne", color="seurat_clusters")
scv.pl.velocity_embedding_stream(adata, basis="tsne", color="seurat_clusters")


#CPM
scv.pl.velocity(adata, var_names=['Tbx1', 'Hoxb1', 'Aplnr', 'Meox1', 'Prdm1', 'Nrg1', 'Gbx2', 'Kdr', 'Grem1', 'Slit1', 'Sema3e', 'Spsb1', 'Chn2', 'Cldn11', 'Gfra1'])
#CPM/AP
scv.pl.velocity(adata, var_names=['Fgf8', 'Ret', 'Hopx', 'Epha2', 'Tjp1'])
#CPM/VP
scv.pl.velocity(adata, var_names=['Tbx5', 'Osr1', 'Wnt2', 'Snai1', 'Aldh1a2'])
#CPM/Ap+VP
scv.pl.velocity(adata, var_names=['Ntng1', 'Gli1', 'Nkx2-5', 'Bmp4'])
#hand1
scv.pl.velocity(adata, var_names=['Hand1'])


#Identify important genes
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(20)

          0         1         2              3         4        5
0     Nr2f2     Mylk3       Ptn          Tmtc1    Slc8a1      Tek
1     Foxp2       Hck     Epha3         Ctnna2    Sh3bgr     Cdh5
2   Cyp3a13     Naa16     Rbms3         Sema5a     Acta2   Adgrl4
3    Pdgfra       Gyg    Sema3a  D830044D21Rik        Qk     Tph2
4   Trim30a    Sema3e      Ror1           Frzb    Synpo2     Esam
5     Rev3l    Hapln1    Sphkap          Itih5    Sphkap   Il20ra
6     Nrxn1     Smyd1      Dkk2          Pitx2   Tgfb1i1   Adgrf5
7      Ets1    Acot12     Pde3a           Dab1      Myh6      Kdr
8    Pdlim5     Fgf10    Plagl1         Igfbp3      Nebl     Tjp1
9      Pid1  Atp6ap1l  Arhgap29        Fam19a4     Pcdh7     Emcn
10    Moxd1     Armh4     Celf2          Negr1       Dsp    Actg2
11     Fat3      Cd36    Mpped2         Hs6st2  Atp6ap1l  Galnt18
12      Nrk    Unc45b      Amot          Vgll3   Gm11592     Flt4
13    Epha4   Tgfb1i1     Gap43         Slc2a9    Igfbp5    Mmrn2
14     Scg5     Rcsd1     Flrt2            Ddc       Mme    Klhl4
15    Rdh10      Cdh6     Smoc2         Ccdc80  Ppp1r14c     Rhoj
16    Csrp2     Nabp1      Nexn          Cpne4   Ccdc141  Gm10801
17      Fn1      Lifr      Cnr1         Snap91      Amph     Nrp1
18     Myo6    Efemp1      Igf1          Abcg1       Ddc    Myzap
19   Cxcl12    Nkx2-6   Ccdc141          Snai2    Trim55     Car8


#Cell cycle
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])
               
#Speed and coherence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

df = adata.obs.groupby('seurat_clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
print(df)

seurat_clusters               0           1           2           3            4           5
velocity_length      947.435425  731.294006  883.488892  697.034790  1139.433838  508.858734
velocity_confidence    0.903120    0.857498    0.901174    0.853661     0.864939    0.850328


#velocity graph
scv.pl.velocity_graph(adata, threshold=.1)

#the graph can be used to draw descendents/anscestors coming from a specified cell. 
x, y = scv.utils.get_cell_transitions(adata, basis='tsne', starting_cell=50)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

#velocity pseudotime - after inferring a distribution over root cells from the graph, it measures the average number of steps it takes to reach a cell after walking along the graph starting from the root cells.
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')


#Page velocity graph
# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='seurat_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
print(df)

     0     1     2    3     4    5
0  0.0  0.00  0.24  0.0  0.04  0.0
1  0.0  0.00  0.13  0.0  0.00  0.0
2  0.0  0.00  0.00  0.0  0.00  0.0
3  0.0  0.01  0.00  0.0  0.00  0.0
4  0.0  0.00  0.00  0.0  0.00  0.0
5  0.0  0.00  0.00  0.0  0.00  0.0

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5)




