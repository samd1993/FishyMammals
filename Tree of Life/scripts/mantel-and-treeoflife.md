# Song et al. (2020) itol mantel figure reproduce attempt 

## After hours and hours of blood sweat and tears I have finally got python, jupyter, and qiime2 to work in peaceful coexistence

### Some notes to remember:

1. Download qiime2 in ubuntu FIRST then when u activate your environment, download ipython and jupyter stuff and kernel in the virtual env. Then there should be code to add kernel to environment. 
2. Then you launch jupyter from the venv and qiime2 can finally be imported in jupyter. This took me forever to figure out. Also had trouble with other packages, but after doing the !sys pip install method in jupyter I had almost no issues. 
3. TreeTime is what I shouldve been using all along, just google TreeTime database to make the tree and then you can use their code to line up all tree, distance matrices, and metadata. Code is so much quicker and easier than in R/excel


```python
import sys
!{sys.executable} -m pip install ecopy
```

    Requirement already satisfied: ecopy in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (0.1.2.2)
    Requirement already satisfied: matplotlib>=1.3.1 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from ecopy) (3.4.1)
    Requirement already satisfied: pandas>=0.13 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from ecopy) (1.2.4)
    Requirement already satisfied: numpy>=1.7 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from ecopy) (1.20.2)
    Requirement already satisfied: patsy>=0.3.0 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from ecopy) (0.5.1)
    Requirement already satisfied: scipy>=0.14 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from ecopy) (1.6.2)
    Requirement already satisfied: pyparsing>=2.2.1 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from matplotlib>=1.3.1->ecopy) (2.4.7)
    Requirement already satisfied: cycler>=0.10 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from matplotlib>=1.3.1->ecopy) (0.10.0)
    Requirement already satisfied: kiwisolver>=1.0.1 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from matplotlib>=1.3.1->ecopy) (1.3.1)
    Requirement already satisfied: pillow>=6.2.0 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from matplotlib>=1.3.1->ecopy) (8.1.2)
    Requirement already satisfied: python-dateutil>=2.7 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from matplotlib>=1.3.1->ecopy) (2.8.1)
    Requirement already satisfied: six in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from cycler>=0.10->matplotlib>=1.3.1->ecopy) (1.15.0)
    Requirement already satisfied: pytz>=2017.3 in /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages (from pandas>=0.13->ecopy) (2021.1)



```python
import skbio as skb
import numpy as np
import pandas as pd
import seaborn as sns
from ecopy import Mantel
```


```python
cd ../

```

    /mnt/c/Users/samde



```python
from matplotlib import pyplot as plt
import pylab as pl
from qiime2 import Artifact
from skbio import TreeNode
from skbio.stats.distance import mantel
from scipy.stats import linregress
from scipy.spatial.distance import squareform, pdist
from os.path import abspath, join
from os import makedirs
import matplotlib as mpl
```


```python
%matplotlib inline
```


```python
pwd
```




    '/mnt/c/Users/samde'




```python
tree_dir = abspath('../samde/python/trees')
host_tree_fp = join(tree_dir, 'timetree.nwk')
host_tree_fp
```




    '/mnt/c/Users/samde/python/trees/timetree.nwk'




```python
host_tree= skb.io.read(host_tree_fp, format='newick', 
                       into=TreeNode,
                       convert_underscores=False)

host_tips = [x.name for x in host_tree.tips()]
len(host_tips)
```




    199




```python
md_dir = '../samde/python/metadata'
host_md_fp = join(md_dir, 'allmerged-metadata4.txt')
host_md = pd.read_csv(host_md_fp, sep='\t')

host_md = host_md.loc[(host_md['treetax'].isin(host_tips))]
```


```python
len(host_md['sampleid'])

```




    887




```python
# filter to just tetrapods
include_classes = ['Mammalia',
                   'Actinopterygii']

host_md = host_md.loc[host_md['host'].isin(include_classes)]
```


```python

dm_dir = abspath('../samde/python/bdiv')
dm_fp = join(dm_dir, 'unweighted_unifrac_distance_matrixdietfixed400.qza')
dm_art = Artifact.load(dm_fp)
dm = dm_art.view(skb.DistanceMatrix)
```


```python
md_ids = set(host_md['sampleid'])

dm_ids = set(dm.ids)

shared_ids = dm_ids & md_ids
```


```python
len(md_ids)
```




    815




```python
len(dm_ids)
```




    786




```python
len(shared_ids)
```




    590




```python
host_md = host_md.loc[host_md['sampleid'].isin(shared_ids)]
dm = dm.filter(list(shared_ids))
```


```python
host_md.shape
```




    (590, 15)




```python
dm.shape
```




    (590, 590)




```python
host_tree = host_tree.shear(host_md['treetax'])
host_tips = [x.name for x in host_tree.tips()]
```


```python
len(host_tips)
```




    108




```python
subset = False
```


```python
# subset the tree if desired
if subset:
    host_subset = np.random.choice(host_tips, size=100, replace=False)
else:
    host_subset = host_tips

host_tree_subset = host_tree.shear(host_subset)

host_ids_subset = host_md.loc[host_md['treetax'].isin(host_subset), 'sampleid']
```


```python
len(host_ids_subset)
```




    590




```python
one_per_sp = True
```


```python
if one_per_sp:
    host_md = host_md.loc[(host_md['sampleid'].isin(host_ids_subset)) &
                          (host_md['sampleid'].isin(dm.ids)),].groupby('treetax').first()

    host_md =  host_md.loc[(host_md['sampleid'].isin(host_ids_subset)) &
                          (host_md['sampleid'].isin(dm.ids)),].groupby('treetax').first().reset_index()
    host_ids_subset = list(set(host_ids_subset) & set(host_md['sampleid']))
```


```python
len(host_ids_subset)
```




    108




```python
patristic_dm = host_tree_subset.tip_tip_distances()
```


```python
patristic_dm.shape
```




    (108, 108)



So now to make diet matrix I went back into relatedness script in R (or start fresh one) and you just use dcast function from reshape2 to make diet matrix from metadata, which in my case is basically binary (0,1), to make distance matrix


```python
md_dir = '../samde/python/metadata'
host_diet_file = join(md_dir, 'alldietd.csv')
host_diet_df = pd.read_csv(host_diet_file, sep=',',index_col=0)
#host_diet_df = host_diet_df.dropna()
host_diet_df.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>carnivore</th>
      <th>corallivore</th>
      <th>detritivore</th>
      <th>herbivore</th>
      <th>omnivore</th>
      <th>planktivore</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>sampleid</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>100</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>101</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
#host_diet_df.set_index('sampleid', inplace=True)  dont need it
```


```python
host_diet_df = host_diet_df.loc[host_ids_subset,:].dropna()
host_diet_df
#for some reason this works...guess i dont need sampleid as a column name like in github
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>carnivore</th>
      <th>corallivore</th>
      <th>detritivore</th>
      <th>herbivore</th>
      <th>omnivore</th>
      <th>planktivore</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>36</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>102Kulan</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>333GraySeal</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>283Koala</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>412Lion</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>101Horse</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>30BeechMarten</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>D1</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>203RedDeer</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>340WesternLowlandGorilla</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>108 rows × 6 columns</p>
</div>




```python
diet_dm = skb.DistanceMatrix(squareform(pdist(host_diet_df.iloc[:, :], metric='jaccard')))
diet_dm.ids = host_diet_df.index
```


```python
diet_dm.shape
```




    (108, 108)




```python
rename_df = host_md.loc[host_md['sampleid'].isin(diet_dm.ids),
                             ['sampleid','treetax']].set_index('sampleid')
rename = [rename_df.loc[x, 'treetax'] for x in diet_dm.ids]
```


```python
diet_dm.ids = rename
```

## Subset DM to remaining hosts and rename by tree


```python
dm_subset = dm.filter(list(host_ids_subset))
```


```python
dm_subset.shape
```




    (108, 108)




```python
rename_df = host_md.loc[host_md['sampleid'].isin(dm_subset.ids),
                             ['sampleid','treetax']].set_index('sampleid')
rename = [rename_df.loc[x, 'treetax'] for x in dm_subset.ids]
```


```python
dm_subset.ids = rename
```


```python
set(dm_subset.ids) == set(patristic_dm.ids)
```




    True




```python
def recursive_mantel(tree, dm1, dm2, min_size=7, **kwargs):
    i = 0
    nodes = tree.count(tips=True) * 2 - 1
    tree.assign_ids()
    node_dict = {}
    for n in tree.postorder():
        label = n.name
        i += 1
        # see if you're below minimum clade size
        if n.count(tips=True) < min_size:
            continue
        # if not, run mantel
        else:
            tips = [t.name for t in n.tips()]
            dm1_s = dm1.filter(tips)
            dm2_s = dm2.filter(tips)
            corr, p, _ = mantel(dm1_s, dm2_s, **kwargs)
            node_dict[n.name] = (corr, p)
    
        pct = (i / nodes)
        
        if np.round(pct * 100) % 10 == 0:
            print(pct)
    return(node_dict)
```


```python
def recursive_mantel_ecopy(tree, dm1, dm2, min_size=7, **kwargs):
    i = 0
    nodes = tree.count(tips=True) * 2 - 1
    tree.assign_ids()
    node_dict = {}
    for n in tree.postorder():
        label = n.name
        i += 1
        # see if you're below minimum clade size
        if n.count(tips=True) < min_size:
            continue
        # if not, run mantel
        else:
            tips = [t.name for t in n.tips()]
            dm1_s = dm1.filter(tips)
            dm2_s = dm2.filter(tips)
            res = Mantel(dm1_s.data, dm2_s.data, **kwargs)
            
            node_dict[n.name] = (res.r_obs, res.pval)
    
        pct = (i / nodes)
        
        if np.round(pct * 100) % 10 == 0:
            print(pct)
    return(node_dict)
```


```python
def recursive_partial_mantel_ecopy(tree, dm1, dm2, dmc, min_size=7, **kwargs):
    i = 0
    nodes = tree.count(tips=True) * 2 - 1
    tree.assign_ids()
    node_dict = {}
    for n in tree.postorder():
        label = n.name
        i += 1
        # see if you're below minimum clade size
        if n.count(tips=True) < min_size:
            continue
        # if not, run mantel
        else:
            tips = [t.name for t in n.tips()]
            dm1_s = dm1.filter(tips)
            dm2_s = dm2.filter(tips)
            dmc_s = dmc.filter(tips)
            res = Mantel(dm1_s.data, dm2_s.data, d_condition=dmc_s.data, **kwargs)
            
            node_dict[n.name] = (res.r_obs, res.pval)
    
        pct = (i / nodes)
        
        if np.round(pct * 100) % 10 == 0:
            print(pct)
    return(node_dict)
```


```python
# assign internal node names, non-integer

host_tree_subset.assign_ids()

for n in host_tree_subset.postorder():
    if n.is_tip():
        continue
    elif n.name[0].isdigit(): 
        n.name = "node%s" % n.name
```


```python
tips = set(patristic_dm.ids) & set(dm_subset.ids) & set(diet_dm.ids)
```


```python
len(tips)
```




    108




```python
run = True
```


```python
if run:
    node_dict_spearman = recursive_mantel(host_tree_subset, patristic_dm, dm_subset, method='spearman')
```

    0.4046511627906977
    0.6976744186046512
    0.7023255813953488
    0.9953488372093023
    1.0



```python
if run:
    node_dict_pearson = recursive_mantel(host_tree_subset, patristic_dm, dm_subset, method='pearson')
```

    0.4046511627906977
    0.6976744186046512
    0.7023255813953488
    0.9953488372093023
    1.0


## Run Partial Mantels tests conditioned on diet


```python
if run:
    node_dict_spearman_partial = recursive_partial_mantel_ecopy(host_tree_subset.shear(tips),
                                      patristic_dm.filter(tips),
                                      dm_subset.filter(tips),
                                      diet_dm.filter(tips),
                                      test='spearman')
```

    /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages/ecopy/matrix_comp/mantel.py:193: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
    To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
      params = np.linalg.lstsq(Z, y_flat)[0]


    0.4046511627906977


    /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages/ecopy/matrix_comp/mantel.py:146: RuntimeWarning: invalid value encountered in true_divide
      y_flat = (y_flat - y_flat.mean())/y_flat.std(ddof=1)
    /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages/ecopy/matrix_comp/mantel.py:145: RuntimeWarning: invalid value encountered in true_divide
      x_flat = (x_flat - x_flat.mean())/x_flat.std(ddof=1)


    0.6976744186046512
    0.7023255813953488
    0.9953488372093023
    1.0



```python
if run:
    node_dict_pearson_partial = recursive_partial_mantel_ecopy(host_tree_subset.shear(tips),
                                      patristic_dm.filter(tips),
                                      dm_subset.filter(tips),
                                      diet_dm.filter(tips),
                                      test='pearson')
```

    0.4046511627906977
    0.6976744186046512
    0.7023255813953488
    0.9953488372093023
    1.0



```python
pwd
```




    '/mnt/c/Users/samde'




```python
out_dir = '../samde/python/mantel'
makedirs(out_dir, exist_ok=True)
out_dir
```




    '../samde/python/mantel'




```python
if run:
    host_tree_subset.write(join(out_dir,'annotated_host_tree.tre'), format='newick')
```


```python
if run:
    node_df_pearson = pd.DataFrame.from_dict(node_dict_pearson, orient='index')
    node_df_pearson.columns = ['r','p']
    node_df_pearson.to_csv(join(out_dir, 'mantel.uwdist.pearson.csv'))
```


```python
if run:
    node_df_spearman =  pd.DataFrame.from_dict(node_dict_spearman, orient='index')
    node_df_spearman.columns = ['r','p']
    node_df_spearman.to_csv(join(out_dir, 'mantel.uwdist.spearman.csv'))
```


```python
if run:
    node_df_pearson_partial = pd.DataFrame.from_dict(node_dict_pearson_partial, orient='index')
    node_df_pearson_partial.columns = ['r','p']
    node_df_pearson_partial.to_csv(join(out_dir, 'mantel_partial.juwdist.pearson.csv'))
```


```python
if run:
    node_df_spearman_partial =  pd.DataFrame.from_dict(node_dict_spearman_partial, orient='index')
    node_df_spearman_partial.columns = ['r','p']
    node_df_spearman_partial.to_csv(join(out_dir, 'mantel_partial.uwdist.spearman.csv'))
```


```python
host_tree_subset = TreeNode.read(join(out_dir,'annotated_host_tree.tre'))
```


```python
node_df_pearson = pd.read_csv(join(out_dir, 'mantel.uwdist.pearson.csv'), index_col=0)
node_df_spearman = pd.read_csv(join(out_dir, 'mantel.uwdist.spearman.csv'), index_col=0)
node_df_pearson_partial = pd.read_csv(join(out_dir, 'mantel_partial.uwdist.pearson.csv'), index_col=0)
node_df_spearman_partial = pd.read_csv(join(out_dir, 'mantel_partial.uwdist.spearman.csv'), index_col=0)
```


```python
def node_values(tips, tree, node_df, tree_dm, dist_dm):
    node_name = tree.lca(tips).name
    
    r = node_df_pearson.loc[node_name, 'r']
    p = node_df_pearson.loc[node_name, 'p']
    
    pats = patristic_dm.filter(tips)
    dsts = dm_subset.filter(tips)

    x = pats.to_series()
    y = dsts.to_series()

    slope, _, _, _, _ = linregress(x, y)
    
    plot = sns.lmplot('x',
                      'y',
                      pd.DataFrame({'x': x, 'y': y}).reset_index(),
                      scatter_kws={'alpha':0.1})
    plt.ylim(0.5, 1.05)

    print(("r: {r}\n"
       "p: {p}\n"
       "m: {m:8.7f}\n").format(r=r, p=p, m=slope))
    
#     plot
```


```python
tips = set(host_md.loc[host_md['host'] == 'Mammalia','treetax']) & set(host_subset)

node_values(tips, host_tree_subset, node_df_pearson, patristic_dm, dm_subset)
```

    /home/samd1993/miniconda3/envs/qiime2-2021.4/lib/python3.8/site-packages/seaborn/_decorators.py:36: FutureWarning: Pass the following variables as keyword args: x, y, data. From version 0.12, the only valid positional argument will be `data`, and passing other arguments without an explicit keyword will result in an error or misinterpretation.
      warnings.warn(


    r: 0.2627149046303847
    p: 0.001
    m: 0.0003126
    



```python
itol_dir = '../samde/python/itol'
```


```python

import matplotlib as mpl
```


```python
node_df_pearson['r'].max()
```




    0.7986208431028278




```python
node_df_spearman['r'].max()
```




    0.7673862831941181




```python

node_df_pearson_partial['r'].max()
```




    0.7986208431028278




```python
node_df_spearman_partial['r'].max()
```




    0.7403988768089856




```python
# create dummy invisible image
# (use the colormap you want to have on the colorbar)
cmax = node_df_pearson['r'].max()
img = plt.imshow(np.array([[0,cmax]]), cmap=mpl.cm.plasma)
img.set_visible(False)

plt.colorbar(orientation="vertical")

# add any other things you want to the figure.
#plt.plot(np.random.rand(30))
```




    <matplotlib.colorbar.Colorbar at 0x7ff6f1258be0>




```python
def colorbar(fp,
             cmax=1,
             size=(4,.7)):
    a = np.array([[0,cmax]])
    pl.figure(figsize=size)
    img = pl.imshow(a, cmap=mpl.cm.plasma)
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.4, 0.8, 0.5])
    pl.colorbar(orientation='horizontal', cax=cax)
    pl.savefig(fp)
```


```python
colorbar(join(itol_dir, 'mantel.uwdist.pearson.itol-r_colorbar.pdf'),
         cmax=node_df_pearson['r'].max(),
         size=(2.5,.6))
```


```python
colorbar(join(itol_dir, 'mantel.uwdist.spearman.itol-r_colorbar.pdf'),
         cmax=node_df_spearman['r'].max(),
         size=(2.5,.6))
```


```python
colorbar(join(itol_dir, 'mantel_partial.uwdist.pearson.itol-r_colorbar.pdf'),
         cmax=node_df_pearson_partial['r'].max(),
         size=(2.5,.6))
```


```python
colorbar(join(itol_dir, 'mantel_partial.uwdist.spearman.itol-r_colorbar.pdf'),
         cmax=node_df_spearman_partial['r'].max(),
         size=(2.5,.6))
```


```python
def get_colors(vals, col, cmap=mpl.cm.plasma):
    norm = mpl.colors.Normalize(vmin=0, vmax=np.max(vals))
    
    colors = []
    
    for i in vals:
        rgb = mpl.cm.ScalarMappable(norm=norm,
                                    cmap=cmap).to_rgba(i)
        hex_val = mpl.colors.rgb2hex(rgb)
        colors.append(hex_val)

    return(colors)

def write_itol_colors(node_df, col, fp, label='branch_colors'):
    out = ('DATASET_STYLE\n\n'
           'SEPARATOR SPACE\n\n'
           'COLOR #ffff00\n\n'
           'DATASET_LABEL %s\n\n'
           'DATA\n'
           '#NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR\n' % label)

    for i, row in node_df.iterrows():
        out += '{0} branch node {1} 2 normal\n'.format(i, row[col])

    with open(fp, 'w') as f:
        f.write(out)
```


```python
#LEGEND_TITLE,Dataset legend
#LEGEND_SHAPES,1,2,3
#LEGEND_COLORS,#ff0000,#00ff00,#0000ff
#LEGEND_LABELS,value1,value2,value3
```


```python
node_df_pearson['r_color'] = get_colors(node_df_pearson['r'], col='r')
node_df_spearman['r_color'] = get_colors(node_df_spearman['r'], col='r')
```


```python
node_df_pearson_partial['r_color'] = get_colors(node_df_pearson_partial['r'], col='r')
node_df_spearman_partial['r_color'] = get_colors(node_df_spearman_partial['r'], col='r')
```


```python
write_itol_colors(node_df_pearson, 
                  'r_color', 
                  join(itol_dir, 'mantel.uwdist.pearson.itol-r_colors.txt'),
                  label='pearson_r')
write_itol_colors(node_df_spearman, 
                  'r_color',
                  join(itol_dir, 'mantel.uwdist.spearman.itol-r_colors.txt'),
                  label='spearman_r')
```


```python
write_itol_colors(node_df_pearson_partial, 
                  'r_color', 
                  join(itol_dir, 'mantel_partial.uwdist.pearson.itol-r_colors.txt'),
                  label='pearson_r_partial')
write_itol_colors(node_df_spearman_partial, 
                  'r_color',
                  join(itol_dir, 'mantel_partial.uwdist.spearman.itol-r_colors.txt'),
                  label='spearman_r_partial')
```

Write pie chart p val indicator file


```python

# 9132,0,50,11000
def write_pval_pies(node_df, col, fp, label='p-vals', breaks = [0.001, 0.01, 0.05]):
    header = ("DATASET_PIECHART\n"
              "SEPARATOR SPACE\n"
              "DATASET_LABEL %s\n"
              "COLOR #ff0000\n"
              "FIELD_COLORS #ff0000\n"
              "FIELD_LABELS p-value\n"
              "MAXIMUM_SIZE 30\n"
              "BORDER_WIDTH 0\n"
              "DATA\n" % label)
    
    break_sizes = [int(30 - i/len(breaks)*30) for i, x in enumerate(breaks)]
    pies = []
    for i, row in node_df.iterrows():
        size = 0
        for s, b in enumerate(breaks):
            if row[col] > b:
                continue
            else:
                size = break_sizes[s]
                break
        if size > 0:
            pies.append("{0} 0 {1} 1".format(i, size))
    
    out = header + '\n'.join(pies)
    
    with open(fp, 'w') as f:
        f.write(out)
```


```python

write_pval_pies(node_df_pearson, 
                'p', 
                join(itol_dir, 'mantel.uwdist.pearson.itol-p_pies.txt'), 
                label='pearson_p')
write_pval_pies(node_df_spearman, 
                'p', 
                join(itol_dir, 'mantel.uwdist.spearman.itol-p_pies.txt'),
                label='spearman_p')
```


```python
write_pval_pies(node_df_pearson_partial, 
                'p',
                join(itol_dir, 'mantel_partial.uwdist.pearson.itol-p_pies.txt'),
                label='pearson_p_partial')
write_pval_pies(node_df_spearman_partial, 
                'p',
                join(itol_dir, 'mantel_partial.uwdist.spearman.itol-p_pies.txt'),
                label='spearman_p_partial')
```


```python
def write_clade_highlights(md_df, tree, tax_col, tree_col, fp,
                           label='clade_highlights',
                           colors=None, alpha=0.4, cmap=mpl.cm.Set3):
    out = ('DATASET_COLORSTRIP\n\n'
           'SEPARATOR SPACE\n\n'
           'COLOR #b2df8a\n\n'
           'DATASET_LABEL %s\n\n'% label)
    
    taxa = md_df[tax_col].unique()
    
    if colors is None:
        clist = cmap(range(len(taxa)), alpha=alpha)
        rgblist = ['rgba(%s,%s,%s,%s)' % (int(x[0]*255),
                                          int(x[1]*255),
                                          int(x[2]*255),
                                          x[3]) for x in clist]
        colors = {x: rgblist[i] for i, x in enumerate(taxa)}
    
    legend_shapes = ' '.join(['1']*len(taxa))
    legend_colors = ''
    legend_labels = ''
    
    for x in colors:
        legend_colors += ' %s' % colors[x]
        legend_labels += ' %s' % x
    
    out += ('LEGEND_TITLE host_colors\n'
            'LEGEND_SHAPES {0}\n'
            'LEGEND_COLORS{1}\n'
            'LEGEND_LABELS{2}\n\nDATA\n').format(legend_shapes,
                                            legend_colors,
                                            legend_labels)
        
    for taxon in taxa:
        tips = set(md_df.loc[md_df[tax_col] == taxon, tree_col])
        node = tree.lowest_common_ancestor(tips)
        
        out += "{0} {1}\n".format(node.name, colors[taxon])
#         out += "{0} range {1} {2}\n".format(node.name, colors[taxon], taxon)
        
    with open(fp, 'w') as f:
        f.write(out)
```


```python
class_colors = {'Actinopterygii': '#b2df8a',
                 'Mammalia': '#1f78b4',}

write_clade_highlights(host_md, 
                       host_tree_subset, 
                       'host', 
                       'treetax',
                       join(itol_dir, 'mantel.itol-Class_colors.txt'),
                       colors=class_colors)
```

## Taxa bar chart for itol


```python
md_dir = '../samde/python/metadata'
tax_md_fp = join(md_dir, 'taxabarchart.csv')
tax_itol = pd.read_csv(tax_md_fp, sep=',')
tax_itol
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>taxtree</th>
      <th>Proteobacteria</th>
      <th>Planctomycetota</th>
      <th>Verrucomicrobiota</th>
      <th>Firmicutes</th>
      <th>Actinobacteriota</th>
      <th>Desulfobacterota</th>
      <th>Bacteroidota</th>
      <th>Cyanobacteria</th>
      <th>Fusobacteriota</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Antidorcas_marsupialis</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>145</td>
      <td>1</td>
      <td>0</td>
      <td>40</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Bos_javanicus</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>70</td>
      <td>0</td>
      <td>0</td>
      <td>21</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Budorcas_taxicolor</td>
      <td>1</td>
      <td>0</td>
      <td>4</td>
      <td>78</td>
      <td>0</td>
      <td>1</td>
      <td>70</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gazella_spekei</td>
      <td>1</td>
      <td>0</td>
      <td>6</td>
      <td>131</td>
      <td>2</td>
      <td>0</td>
      <td>32</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Ovis_ammon</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>194</td>
      <td>306</td>
      <td>0</td>
      <td>34</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>195</th>
      <td>Ovis_aries</td>
      <td>374</td>
      <td>22</td>
      <td>513</td>
      <td>56875</td>
      <td>167</td>
      <td>233</td>
      <td>23001</td>
      <td>608</td>
      <td>0</td>
    </tr>
    <tr>
      <th>196</th>
      <td>Porites_lobata</td>
      <td>186009</td>
      <td>237</td>
      <td>175</td>
      <td>51</td>
      <td>36</td>
      <td>48</td>
      <td>2360</td>
      <td>685</td>
      <td>93</td>
    </tr>
    <tr>
      <th>197</th>
      <td>Pocillopora_elegans</td>
      <td>183608</td>
      <td>293</td>
      <td>228</td>
      <td>21</td>
      <td>90</td>
      <td>149</td>
      <td>3019</td>
      <td>9627</td>
      <td>168</td>
    </tr>
    <tr>
      <th>198</th>
      <td>Acropora_hyacinthus</td>
      <td>39601</td>
      <td>817</td>
      <td>1547</td>
      <td>33</td>
      <td>129</td>
      <td>140</td>
      <td>5694</td>
      <td>3405</td>
      <td>91</td>
    </tr>
    <tr>
      <th>199</th>
      <td>Stegastes_nigricans</td>
      <td>134862</td>
      <td>47704</td>
      <td>16916</td>
      <td>22092</td>
      <td>6293</td>
      <td>7938</td>
      <td>24874</td>
      <td>13554</td>
      <td>12723</td>
    </tr>
  </tbody>
</table>
<p>200 rows × 10 columns</p>
</div>




```python
def get_qualitative_colors(fields, cmap=mpl.cm.Set3, alpha=0.4):
    clist = cmap(range(len(taxa)), alpha=alpha)
    rgblist = ['rgba(%s,%s,%s,%s)' % (int(x[0]*255),
                                      int(x[1]*255),
                                      int(x[2]*255),
                                      x[3]) for x in clist]
    return(rgblist)

def format_iTOL_multibar(fields, md, 
                         tree_ref_col=None, 
                         field_colors=None, 
                         field_labels=None, 
                         dataset_label='Multibar Chart', 
                         dataset_color=None, 
                         legend=True, 
                         unstacked=False,
                         width=1000, 
                         margin=0,
                         alpha=0.8,
                        ):
    """
    fields: array of columns titles in [metadata] to chart
    metadata: pd.df containing samples to graph and data
    """
    
    if field_labels is None:
        field_labels=fields
    
    if field_colors is None:
        field_colors=get_qualitative_colors(fields, alpha=alpha)
    
    if tree_ref_col is None:
        tree_ref_col=metadata.columns[0]

    if dataset_color is None:
        dataset_color = "#00FF00"

    # remove nans
    metadata = md.loc[:,fields + [tree_ref_col]].dropna()
    
    outstring = ''
    
    outstring += 'DATASET_MULTIBAR\n'
    outstring += 'SEPARATOR TAB\n'
    outstring += 'DATASET_LABEL\t%s\n' % dataset_label
    outstring += 'COLOR\t%s\n' % dataset_color
    outstring += 'FIELD_COLORS\t%s\n' % '\t'.join(field_colors)
    outstring += 'FIELD_LABELS\t%s\n' % '\t'.join(field_labels)    
    
    if legend:
        outstring += 'LEGEND_TITLE\tDataset legend\n'
        outstring += 'LEGEND_SHAPES\t%s\n' % '\t'.join(['1']*len(fields))
        outstring += 'LEGEND_COLORS\t%s\n' % '\t'.join(field_colors)
        outstring += 'LEGEND_LABELS\t%s\n' % '\t'.join(field_labels)
    
    outstring += 'MARGIN\t%s\n' % margin
    outstring += 'WIDTH\t%s\n' % width

    if unstacked:
        outstring += 'ALIGN_FIELDS\t1\n'
        
    outstring += 'DATA\n'
    for index, row in metadata.iterrows():
        outstring += row[tree_ref_col].replace(' ', '_') + '\t%s\n' % '\t'.join([str(row[x]) for x in fields])

    return(outstring)
```


```python
fields = ['Proteobacteria',
 'Planctomycetota',
 'Verrucomicrobiota',
 'Firmicutes',
 'Actinobacteriota',
 'Desulfobacterota',
 'Bacteroidota',
 'Cyanobacteria',
 'Fusobacteriota']

```


```python
bar = format_iTOL_multibar(fields, 
                     tax_itol, 
                     tree_ref_col='taxtree', 
                     field_colors=['#a6cee3',
                                   '#cab2d6',
                                   '#1f78b4',
                                   '#33a02c',
                                   '#6a3d9a',
                                   '#b2df8a',
                                   '#fb9a99',
                                   '#e31a1c',
                                   '#ff7f00',], 
                     field_labels=['Proteobacteria',
 'Planctomycetota',
 'Verrucomicrobiota',
 'Firmicutes',
 'Actinobacteriota',
 'Desulfobacterota',
 'Bacteroidota',
 'Cyanobacteria',
 'Fusobacteriota'], 
                     dataset_label='Taxa barchart2', 
                     dataset_color=None, 
                     legend=True,
                     width=200, 
                     alpha=1.0)

with open(join(itol_dir, 'mantel.itol-taxa_bar.txt'), 'w') as f:
    f.write(bar)
```

Now to make the diet highlights for tree


```python
def write_diet_highlights(md_df, tree, tax_col, tree_col, fp,
                           label='diet_colors',
                           colors=None, alpha=0.4, cmap=mpl.cm.Set3):
    out = ('DATASET_COLORSTRIP\n\n'
           'SEPARATOR SPACE\n\n'
           'COLOR #b2df8a\n\n'
           'DATASET_LABEL %s\n\n'% label)
    
    taxa = md_df[tax_col].unique()
    
    if colors is None:
        clist = cmap(range(len(taxa)), alpha=alpha)
        rgblist = ['rgba(%s,%s,%s,%s)' % (int(x[0]*255),
                                          int(x[1]*255),
                                          int(x[2]*255),
                                          x[3]) for x in clist]
        colors = {x: rgblist[i] for i, x in enumerate(taxa)}
    
    legend_shapes = ' '.join(['1']*len(taxa))
    legend_colors = ''
    legend_labels = ''
    
    for x in colors:
        legend_colors += ' %s' % colors[x]
        legend_labels += ' %s' % x
    
    out += ('LEGEND_TITLE diet_colors\n'
            'LEGEND_SHAPES {0}\n'
            'LEGEND_COLORS{1}\n'
            'LEGEND_LABELS{2}\n\nDATA\n').format(legend_shapes,
                                            legend_colors,
                                            legend_labels)
        
    for taxon in taxa:
        tips = set(md_df.loc[md_df[tax_col] == taxon, tree_col])
        node = tree.lowest_common_ancestor(tips)
        
        out += "{0} {1}\n".format(node.name, colors[taxon])
#         out += "{0} range {1} {2}\n".format(node.name, colors[taxon], taxon)
        
    with open(fp, 'w') as f:
        f.write(out)
```


```python
diet_colors = {'carnivore': '#b2df8a',
                 'herbivore': '#1f78b4',
                   'omnivore' : '#33a02c',
                      'planktivore' : '#6a3d9a',
                        'detritivore' : '#b2df8a',
                          'corallivore' : '#fb9a99'}

write_diet_highlights(host_md, 
                       host_tree_subset, 
                       'diet2', 
                       'treetax',
                       join(itol_dir, 'mantel.itol-diet_colors.txt'),
                       colors=diet_colors)
```


```python
host_md
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>treetax</th>
      <th>sampleid</th>
      <th>sampletaxid</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>diet</th>
      <th>run</th>
      <th>diet1</th>
      <th>host</th>
      <th>diet2</th>
      <th>diet3</th>
      <th>type</th>
      <th>taxaname</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Abudefduf_sexfasciatus</td>
      <td>25</td>
      <td>Abudefduf_sexfasciatus</td>
      <td>Ovalentaria</td>
      <td>Pomacentridae</td>
      <td>Abudefduf</td>
      <td>sexfasciatus</td>
      <td>planktivore</td>
      <td>c</td>
      <td>omnivore</td>
      <td>Actinopterygii</td>
      <td>planktivore</td>
      <td>planktivore</td>
      <td>animal</td>
      <td>"Abudefduf sexfasciatus"</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Aepyceros_melampus</td>
      <td>410Impala</td>
      <td>Aepyceros_melampus</td>
      <td>Artiodactyla</td>
      <td>Bovidae</td>
      <td>Aepyceros</td>
      <td>melampus</td>
      <td>H</td>
      <td>d</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Aepyceros melampus"</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Ailuropoda_melanoleuca</td>
      <td>GP</td>
      <td>Ailuropoda_melanoleuca</td>
      <td>Carnivora</td>
      <td>Ursidae</td>
      <td>Ailuropoda</td>
      <td>melanoleuca</td>
      <td>H</td>
      <td>a</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Ailuropoda melanoleuca"</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Ailurus_fulgens</td>
      <td>RP</td>
      <td>Ailurus_fulgens</td>
      <td>Carnivora</td>
      <td>Ursidae</td>
      <td>Ailurus</td>
      <td>fulgens</td>
      <td>H</td>
      <td>a</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Ailurus fulgens"</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Apodemus_flavicollis</td>
      <td>21YellowneckedFieldMouse</td>
      <td>Apodemus_flavicollis</td>
      <td>Rodentia</td>
      <td>Muridae</td>
      <td>Apodemus</td>
      <td>flavicollis</td>
      <td>O</td>
      <td>d</td>
      <td>omnivore</td>
      <td>Mammalia</td>
      <td>omnivore</td>
      <td>O</td>
      <td>animal</td>
      <td>"Apodemus flavicollis"</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>103</th>
      <td>Trachypithecus_hatinhensis</td>
      <td>238HanumanLangur</td>
      <td>Trachypithecus_hatinhensis</td>
      <td>Primates</td>
      <td>Cercopithecidae</td>
      <td>Trachypithecus</td>
      <td>hatinhensis</td>
      <td>H</td>
      <td>d</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Trachypithecus hatinhensis"</td>
    </tr>
    <tr>
      <th>104</th>
      <td>Trichosurus_vulpecula</td>
      <td>287CommonBrushtail</td>
      <td>Trichosurus_vulpecula</td>
      <td>Diprotodontia</td>
      <td>Phalangeridae</td>
      <td>Trichosurus</td>
      <td>vulpecula</td>
      <td>O</td>
      <td>d</td>
      <td>omnivore</td>
      <td>Mammalia</td>
      <td>omnivore</td>
      <td>O</td>
      <td>animal</td>
      <td>"Trichosurus vulpecula"</td>
    </tr>
    <tr>
      <th>105</th>
      <td>Ursus_arctos</td>
      <td>231BrownBear</td>
      <td>Ursus_arctos</td>
      <td>Carnivora</td>
      <td>Ursidae</td>
      <td>Ursus</td>
      <td>arctos</td>
      <td>O</td>
      <td>d</td>
      <td>omnivore</td>
      <td>Mammalia</td>
      <td>omnivore</td>
      <td>O</td>
      <td>animal</td>
      <td>"Ursus arctos"</td>
    </tr>
    <tr>
      <th>106</th>
      <td>Vulpes_vulpes</td>
      <td>109RedFox</td>
      <td>Vulpes_vulpes</td>
      <td>Carnivora</td>
      <td>Canidae</td>
      <td>Vulpes</td>
      <td>vulpes</td>
      <td>O</td>
      <td>d</td>
      <td>omnivore</td>
      <td>Mammalia</td>
      <td>omnivore</td>
      <td>O</td>
      <td>animal</td>
      <td>"Vulpes vulpes"</td>
    </tr>
    <tr>
      <th>107</th>
      <td>Zebrasoma_scopas</td>
      <td>1</td>
      <td>Zebrasoma_scopas</td>
      <td>Acanthuriformes</td>
      <td>Acanthuridae</td>
      <td>Zebrasoma</td>
      <td>scopas</td>
      <td>herbivore</td>
      <td>c</td>
      <td>herbivore</td>
      <td>Actinopterygii</td>
      <td>herbivore</td>
      <td>herbivore</td>
      <td>animal</td>
      <td>"Zebrasoma scopas"</td>
    </tr>
  </tbody>
</table>
<p>108 rows × 15 columns</p>
</div>




```python
host_tree_subset
```




    <TreeNode, name: node419, internal node count: 106, tips count: 108>



## Now trying to make figure with seaborn


```python
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
```


```python
cd 'OneDrive - UCLA IT Services/Fish Project/ISLAND COMPILED PROJECT/python'
```

    /mnt/c/Users/samde/OneDrive - UCLA IT Services/Fish Project/ISLAND COMPILED PROJECT/python



```python
df = pd.read_csv('faith_pd_merged.txt', delimiter = "\t")
```


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>sampletaxid</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>diet</th>
      <th>run</th>
      <th>diet1</th>
      <th>host</th>
      <th>diet2</th>
      <th>diet3</th>
      <th>type</th>
      <th>taxaname</th>
      <th>tax</th>
      <th>faith_pd</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GP</td>
      <td>Ailuropoda_melanoleuca</td>
      <td>Carnivora</td>
      <td>Ursidae</td>
      <td>Ailuropoda</td>
      <td>melanoleuca</td>
      <td>H</td>
      <td>a</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Ailuropoda melanoleuca"</td>
      <td>Ailuropoda_melanoleuca</td>
      <td>4.394684</td>
    </tr>
    <tr>
      <th>1</th>
      <td>RP</td>
      <td>Ailurus_fulgens</td>
      <td>Carnivora</td>
      <td>Ursidae</td>
      <td>Ailurus</td>
      <td>fulgens</td>
      <td>H</td>
      <td>a</td>
      <td>herbivore</td>
      <td>Mammalia</td>
      <td>herbivore</td>
      <td>H</td>
      <td>animal</td>
      <td>"Ailurus fulgens"</td>
      <td>Ailurus_fulgens</td>
      <td>4.074280</td>
    </tr>
    <tr>
      <th>2</th>
      <td>29</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>Tetraodontiformes</td>
      <td>Balistidae</td>
      <td>Rhinecanthus</td>
      <td>aculeatus</td>
      <td>omnivore</td>
      <td>c</td>
      <td>omnivore</td>
      <td>Actinopterygii</td>
      <td>omnivore</td>
      <td>omnivore</td>
      <td>animal</td>
      <td>"Rhinecanthus aculeatus"</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>25.982524</td>
    </tr>
    <tr>
      <th>3</th>
      <td>48</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>Tetraodontiformes</td>
      <td>Balistidae</td>
      <td>Rhinecanthus</td>
      <td>aculeatus</td>
      <td>omnivore</td>
      <td>c</td>
      <td>omnivore</td>
      <td>Actinopterygii</td>
      <td>omnivore</td>
      <td>omnivore</td>
      <td>animal</td>
      <td>"Rhinecanthus aculeatus"</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>17.935094</td>
    </tr>
    <tr>
      <th>4</th>
      <td>74</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>Tetraodontiformes</td>
      <td>Balistidae</td>
      <td>Rhinecanthus</td>
      <td>aculeatus</td>
      <td>omnivore</td>
      <td>c</td>
      <td>omnivore</td>
      <td>Actinopterygii</td>
      <td>omnivore</td>
      <td>omnivore</td>
      <td>animal</td>
      <td>"Rhinecanthus aculeatus"</td>
      <td>Rhinecanthus_aculeatus</td>
      <td>23.936364</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>781</th>
      <td>M158</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>22.890753</td>
    </tr>
    <tr>
      <th>782</th>
      <td>M251</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>24.949555</td>
    </tr>
    <tr>
      <th>783</th>
      <td>M66</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>22.263455</td>
    </tr>
    <tr>
      <th>784</th>
      <td>M96</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>18.289411</td>
    </tr>
    <tr>
      <th>785</th>
      <td>d</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>6.671163</td>
    </tr>
  </tbody>
</table>
<p>786 rows × 16 columns</p>
</div>




```python
df2 = df[(df["host"] == "Mammalia") | (df["host"] == "Actinopterygii")]
```


```python
df2 = df2[df2['host'].notnull()].copy()
df2['host'] = df2['host'].astype('category')
df2['host']
```




    0            Mammalia
    1            Mammalia
    2      Actinopterygii
    3      Actinopterygii
    4      Actinopterygii
                ...      
    773             algae
    774             algae
    775             algae
    776             algae
    777             algae
    Name: host, Length: 732, dtype: category
    Categories (5, object): ['Actinopterygii', 'Mammalia', 'algae', 'coral', 'water']




```python
#df2['host'].cat.reorder_categories(['Mammalia','Actinopterygii','coral','water','algae'],inplace=True)
```


```python
sns.set(font_scale = 2)
sns.set_style(style='white')
s1=sns.catplot(x="diet1", y="faith_pd", hue="host", kind="box", data=df2,showfliers=False,aspect=4/3,linewidth=3,width=0.5,palette=["#0072B2", "#D55E00"])
s1.set(xlabel=None)

s1.set(ylabel="Alpha Diversity")

plt.show()
```


    
![png](output_111_0.png)
    



```python

```


```python

```


```python
sns.set(font_scale = 1.4)
sns.set_style(style='white')
s2=sns.catplot(x="host", y="faith_pd", kind="box", data=df,showfliers=False,aspect=7/3,linewidth=3,palette=["#0072B2", "#D55E00",
                "gray","gray","gray","gray","gray","gray"])
s2.set(xlabel=None)
s2.set(ylabel="Alpha Diversity")


plt.show()
```


    
![png](output_114_0.png)
    



```python
sns.set(font_scale = 2)
sns.set_style(style='white')
s3=sns.catplot(x="diet1", y="faith_pd", hue="host", kind="box", data=df,showfliers=False,aspect=4/2,linewidth=3,width=0.5)
s3.set(xlabel=None)

s3.set(ylabel="Alpha Diversity")

plt.show()
```


    
![png](output_115_0.png)
    



```python
f, axs = plt.subplots(1,2,
                      figsize=(17,7),
                      sharey=True,gridspec_kw=dict(width_ratios=[3,1.7]))



plt.xticks(rotation=45)
g1=sns.boxplot(x="host", y="faith_pd", data=df,showfliers=False,linewidth=3,palette=["#0072B2", "#D55E00",
                "gray","gray","gray","gray","gray","gray"],ax=axs[0])
g1.set(xlabel=None)
g1.set_xticklabels(g1.get_xticklabels(), rotation=50)
g1.spines['right'].set_visible(False)
g1.spines['top'].set_visible(False)
g1.set_ylabel('Alpha Diversity',size=35,labelpad=20)


g2=sns.boxplot(x="diet1", y="faith_pd", hue="host", data=df2,showfliers=False,linewidth=3,width=0.5,palette=["#0072B2", "#D55E00"],ax=axs[1])
g2.set(xlabel=None)
g2.set(ylabel=None)
g2.set_xticklabels(g2.get_xticklabels(), rotation=50)
g2.spines['right'].set_visible(False)
g2.spines['top'].set_visible(False)
g2.spines['left'].set_visible(False)


sns.set_style(style='white')
plt.tight_layout()
plt.legend(bbox_to_anchor=(0.6, 1), loc='upper left', borderaxespad=0, handlelength=1.5)
sns.set(font_scale = 2)
sns.set_style(style='white')


plt.savefig('AlphaDiversityHostDiet.pdf', dpi=300, bbox_inches='tight')
plt.show()
```


    
![png](output_116_0.png)
    

