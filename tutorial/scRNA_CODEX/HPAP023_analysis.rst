.. code:: ipython3

    import numpy as np
    import pandas as pd
    
    import matplotlib.pyplot as plt
    plt.rcParams["figure.figsize"] = (6, 4)
    
    import anndata as ad
    import scanpy as sc
    import cellink
    from cellink import Cellink

.. code:: ipython3

    # read in rna data
    rna = pd.read_csv('/Users/luosanj/Desktop/project2/data/HPAP023/HPAP-023_exp_mat.csv')

.. code:: ipython3

    rna_feature = rna.iloc[:,0]
    rna_obs = rna.columns[1:]
    rna_mat = rna.iloc[:, 1:].to_numpy().T
    # convert to anndata
    rna_adata = ad.AnnData(rna_mat, dtype = np.float32)
    rna_adata.var_names = rna_feature
    rna_adata.obs_names = rna_obs
    # read in cell type labels of scRNA
    rna_ct = pd.read_csv('/Users/luosanj/Desktop/project2/data/HPAP_cell_annotations/HPAP-023_cell_annotations.csv')
    label_rna = rna_ct['Celltype'].to_numpy()
    rna_adata.obs['cell_type'] = label_rna


.. parsed-literal::

    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/anndata/_core/anndata.py:522: FutureWarning: The dtype argument is deprecated and will be removed in late 2024.
      warnings.warn(


.. code:: ipython3

    protein_adata = sc.read_h5ad('/Users/luosanj/Desktop/project2/data/HPAP023/HPAP023-T1D-17y-7y-F-T_3_2_2-OCT.h5ad')

.. code:: ipython3

    protein_adata = protein_adata[protein_adata.obs['cell_type'] != 'Unidentified'].copy()

.. code:: ipython3

    protein_adata_filter = protein_adata[(protein_adata.obsm['spatial'][:,0] > 5000) & 
                                         (protein_adata.obsm['spatial'][:,1] < 5000)]

.. code:: ipython3

    def map_cell_types(cell_type):
        if cell_type in ['Lymphatic Endothelial Cell', 'Vascular Endothelial Cell']:
            return 'Endothelial'
        elif cell_type in ['T Cell', 'Macrophage Cell', 'Monocyte or Dendritic Cell']:
            return 'Immune'
        elif cell_type in ['Gamma Cell', 'Delta Cell']:
            return 'Delta-PP'
        elif cell_type == 'Exxocrine Cell':
            return 'Exocrine'
        elif cell_type == 'Exocrine Cell':
            return 'Exocrine'
        elif cell_type == 'Duct Cell':
            return 'Ductal'
        elif cell_type == 'Alpha Cell':
            return 'Alpha'
        elif cell_type == 'Acinar':
            return 'Exocrine'
        elif cell_type == 'Beta Cell':
            return 'Beta'
        else:
            return cell_type  # Return the original name if it doesn't match the above
    
    # Apply the mapping function to each cell type in the 'cell_type' column
    rna_adata.obs['cell_type'] = rna_adata.obs['cell_type'].apply(map_cell_types)
    protein_adata_filter.obs['cell_type'] = protein_adata_filter.obs['cell_type'].apply(map_cell_types)


.. parsed-literal::

    /var/folders/zz/067t14gd4nj4836tcg_c5prw0000gt/T/ipykernel_63606/1918363788.py:25: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      protein_adata_filter.obs['cell_type'] = protein_adata_filter.obs['cell_type'].apply(map_cell_types)


.. code:: ipython3

    # # process all RNA features
    sc.pp.normalize_total(rna_adata)
    sc.pp.log1p(rna_adata)
    #sc.pp.scale(rna_adata)
    sc.pp.highly_variable_genes(rna_adata, n_top_genes=5000)
    # # only retain highly variable genes
    rna_adata1 = rna_adata[:, rna_adata.var.highly_variable].copy()
    # # plot UMAPs of rna cells based on all active rna markers
    
    sc.pp.neighbors(rna_adata1, n_neighbors=15)
    sc.tl.umap(rna_adata1)
    sc.pl.umap(rna_adata1, color='cell_type')


.. parsed-literal::

    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/scanpy/preprocessing/_highly_variable_genes.py:226: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
      disp_grouped = df.groupby("mean_bin")["dispersions"]


.. parsed-literal::

    WARNING: You’re trying to run this on 5000 dimensions of `.X`, if you really want this, set `use_rep='X'`.
             Falling back to preprocessing with `sc.pp.pca` and default params.


.. parsed-literal::

    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/anndata/_core/anndata.py:522: FutureWarning: The dtype argument is deprecated and will be removed in late 2024.
      warnings.warn(
    OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/scanpy/plotting/_tools/scatterplots.py:1234: FutureWarning: The default value of 'ignore' for the `na_action` parameter in pandas.Categorical.map is deprecated and will be changed to 'None' in a future version. Please set na_action to the desired value to avoid seeing this warning
      color_vector = pd.Categorical(values.map(color_map))
    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/scanpy/plotting/_tools/scatterplots.py:394: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
      cax = scatter(



.. image:: HPAP023_analysis_files/HPAP023_analysis_7_3.png


.. code:: ipython3

    correspondence = pd.read_csv('/Users/luosanj/Desktop/project2/data/protein_gene_relationship.csv')
    rna_protein_correspondence = []
    
    for i in range(correspondence.shape[0]):
        curr_protein_name, curr_rna_names = correspondence.iloc[i]
        if curr_protein_name not in protein_adata.var_names:
            continue
        if curr_rna_names.find('Ignore') != -1: # some correspondence ignored eg. protein isoform to one gene
            continue
        curr_rna_names = curr_rna_names.split('/') # eg. one protein to multiple genes
        for r in curr_rna_names:
            if r in rna_adata1.var_names:
                rna_protein_correspondence.append([r, curr_protein_name])
                
    rna_protein_correspondence = np.array(rna_protein_correspondence)

.. code:: ipython3

    # correspondence information
    rna_shared = rna_adata1[:, rna_protein_correspondence[:, 0]].copy()
    protein_shared = protein_adata_filter[:, rna_protein_correspondence[:, 1]].copy()


.. parsed-literal::

    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/anndata/_core/anndata.py:1908: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")
    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/anndata/_core/anndata.py:1908: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")


.. code:: ipython3

    extract_ct = ['Alpha', 'Exocrine', 'Delta-PP', 'Ductal']
    rna_shared = rna_shared[rna_shared.obs['cell_type'].isin(extract_ct)]
    protein_shared = protein_shared[protein_shared.obs['cell_type'].isin(extract_ct)]
    combined_adata = rna_adata1[rna_adata1.obs['cell_type'].isin(extract_ct)]
    protein_adata_filter = protein_adata_filter[protein_adata_filter.obs['cell_type'].isin(extract_ct)]

.. code:: ipython3

    # Make sure no column is static
    mask = (
        (rna_shared.X.toarray().std(axis=0) > 0.05) 
        & (protein_shared.X.std(axis=0) > 0.05)
    )
    rna_shared = rna_shared[:, mask].copy()
    protein_shared = protein_shared[:, mask].copy()
    print([rna_shared.shape,protein_shared.shape])


.. parsed-literal::

    [(278, 17), (10642, 17)]


.. parsed-literal::

    /Users/luosanj/opt/anaconda3/envs/alignment/lib/python3.9/site-packages/anndata/_core/anndata.py:1908: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")


.. code:: ipython3

    arr = [rna_shared, protein_shared]
    cellink = Cellink(full_ann1 = combined_adata, full_ann2 = protein_adata_filter, shared_ann1 = rna_shared, shared_ann2 = protein_shared)
    cellink.split_into_batches(arr, 300, seed = 100)


.. parsed-literal::

    Cell annotations are provided. Perform Iteratively OT!
    The first modality is split into 1 batches, and max batch size is 278.
    The second modality is split into 35 batches, and max batch size is 306.
    Batch to batch correspondence is:
      ['0<->0', '0<->1', '0<->2', '0<->3', '0<->4', '0<->5', '0<->6', '0<->7', '0<->8', '0<->9', '0<->10', '0<->11', '0<->12', '0<->13', '0<->14', '0<->15', '0<->16', '0<->17', '0<->18', '0<->19', '0<->20', '0<->21', '0<->22', '0<->23', '0<->24', '0<->25', '0<->26', '0<->27', '0<->28', '0<->29', '0<->30', '0<->31', '0<->32', '0<->33', '0<->34'].


.. code:: ipython3

    cellink.alignment(wt1 = 0.7, wt2 = 0.7, lambd = 0.01, numItermax = 1000, reg = 0.01, 
                      reg_m1 = (20, 0), reg_m2 = (0, 20), iterative = True, sparse = False)


.. parsed-literal::

    Now at batch 0<->0...
    60 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.17%! 
    
    There are 19 unmatched samples and 259 matched samples in data1!
    
    62 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.05%! 
    
    There are 12 unmatched samples and 292 matched samples in data2!
    
    Now at batch 0<->1...
    44 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.73%! 
    
    There are 23 unmatched samples and 255 matched samples in data1!
    
    52 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.71%! 
    
    There are 10 unmatched samples and 294 matched samples in data2!
    
    Now at batch 0<->2...
    41 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    49 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.39%! 
    
    There are 14 unmatched samples and 290 matched samples in data2!
    
    Now at batch 0<->3...
    40 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.01%! 
    
    There are 25 unmatched samples and 253 matched samples in data1!
    
    45 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.07%! 
    
    There are 15 unmatched samples and 289 matched samples in data2!
    
    Now at batch 0<->4...
    45 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.17%! 
    
    There are 19 unmatched samples and 259 matched samples in data1!
    
    56 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.39%! 
    
    There are 14 unmatched samples and 290 matched samples in data2!
    
    Now at batch 0<->5...
    43 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.36999999999999%! 
    
    There are 24 unmatched samples and 254 matched samples in data1!
    
    59 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.42%! 
    
    There are 20 unmatched samples and 284 matched samples in data2!
    
    Now at batch 0<->6...
    52 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.88%! 
    
    There are 17 unmatched samples and 261 matched samples in data1!
    
    51 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 98.36%! 
    
    There are 5 unmatched samples and 299 matched samples in data2!
    
    Now at batch 0<->7...
    48 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.64999999999999%! 
    
    There are 26 unmatched samples and 252 matched samples in data1!
    
    57 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.74000000000001%! 
    
    There are 16 unmatched samples and 288 matched samples in data2!
    
    Now at batch 0<->8...
    45 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.32000000000001%! 
    
    There are 13 unmatched samples and 265 matched samples in data1!
    
    50 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.71%! 
    
    There are 10 unmatched samples and 294 matched samples in data2!
    
    Now at batch 0<->9...
    37 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.17%! 
    
    There are 19 unmatched samples and 259 matched samples in data1!
    
    49 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.07%! 
    
    There are 15 unmatched samples and 289 matched samples in data2!
    
    Now at batch 0<->10...
    59 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.81%! 
    
    There are 20 unmatched samples and 258 matched samples in data1!
    
    55 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.05%! 
    
    There are 12 unmatched samples and 292 matched samples in data2!
    
    Now at batch 0<->11...
    39 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.53%! 
    
    There are 18 unmatched samples and 260 matched samples in data1!
    
    44 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 98.03%! 
    
    There are 6 unmatched samples and 298 matched samples in data2!
    
    Now at batch 0<->12...
    49 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.29%! 
    
    There are 27 unmatched samples and 251 matched samples in data1!
    
    56 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.08%! 
    
    There are 18 unmatched samples and 286 matched samples in data2!
    
    Now at batch 0<->13...
    37 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.29%! 
    
    There are 27 unmatched samples and 251 matched samples in data1!
    
    48 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.07%! 
    
    There are 15 unmatched samples and 289 matched samples in data2!
    
    Now at batch 0<->14...
    69 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.36999999999999%! 
    
    There are 24 unmatched samples and 254 matched samples in data1!
    
    58 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.07%! 
    
    There are 15 unmatched samples and 289 matched samples in data2!
    
    Now at batch 0<->15...
    57 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.45%! 
    
    There are 21 unmatched samples and 257 matched samples in data1!
    
    63 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.07%! 
    
    There are 15 unmatched samples and 289 matched samples in data2!
    
    Now at batch 0<->16...
    48 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.36999999999999%! 
    
    There are 24 unmatched samples and 254 matched samples in data1!
    
    53 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.72%! 
    
    There are 13 unmatched samples and 291 matched samples in data2!
    
    Now at batch 0<->17...
    51 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.73%! 
    
    There are 23 unmatched samples and 255 matched samples in data1!
    
    54 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.71%! 
    
    There are 10 unmatched samples and 294 matched samples in data2!
    
    Now at batch 0<->18...
    46 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.64999999999999%! 
    
    There are 26 unmatched samples and 252 matched samples in data1!
    
    50 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.75999999999999%! 
    
    There are 22 unmatched samples and 282 matched samples in data2!
    
    Now at batch 0<->19...
    50 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    56 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.39%! 
    
    There are 14 unmatched samples and 290 matched samples in data2!
    
    Now at batch 0<->20...
    45 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.64999999999999%! 
    
    There are 26 unmatched samples and 252 matched samples in data1!
    
    53 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.72%! 
    
    There are 13 unmatched samples and 291 matched samples in data2!
    
    Now at batch 0<->21...
    56 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.81%! 
    
    There are 20 unmatched samples and 258 matched samples in data1!
    
    60 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.08%! 
    
    There are 18 unmatched samples and 286 matched samples in data2!
    
    Now at batch 0<->22...
    58 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 93.17%! 
    
    There are 19 unmatched samples and 259 matched samples in data1!
    
    55 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.38%! 
    
    There are 11 unmatched samples and 293 matched samples in data2!
    
    Now at batch 0<->23...
    49 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    49 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 97.04%! 
    
    There are 9 unmatched samples and 295 matched samples in data2!
    
    Now at batch 0<->24...
    47 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.29%! 
    
    There are 27 unmatched samples and 251 matched samples in data1!
    
    60 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.78%! 
    
    There are 25 unmatched samples and 279 matched samples in data2!
    
    Now at batch 0<->25...
    46 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    60 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 96.38%! 
    
    There are 11 unmatched samples and 293 matched samples in data2!
    
    Now at batch 0<->26...
    43 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    53 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.72%! 
    
    There are 13 unmatched samples and 291 matched samples in data2!
    
    Now at batch 0<->27...
    40 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.73%! 
    
    There are 23 unmatched samples and 255 matched samples in data1!
    
    42 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.72%! 
    
    There are 13 unmatched samples and 291 matched samples in data2!
    
    Now at batch 0<->28...
    43 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.09%! 
    
    There are 22 unmatched samples and 256 matched samples in data1!
    
    54 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.74000000000001%! 
    
    There are 16 unmatched samples and 288 matched samples in data2!
    
    Now at batch 0<->29...
    58 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.45%! 
    
    There are 21 unmatched samples and 257 matched samples in data1!
    
    55 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 97.04%! 
    
    There are 9 unmatched samples and 295 matched samples in data2!
    
    Now at batch 0<->30...
    46 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 91.01%! 
    
    There are 25 unmatched samples and 253 matched samples in data1!
    
    52 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.08%! 
    
    There are 18 unmatched samples and 286 matched samples in data2!
    
    Now at batch 0<->31...
    47 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.81%! 
    
    There are 20 unmatched samples and 258 matched samples in data1!
    
    54 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.72%! 
    
    There are 13 unmatched samples and 291 matched samples in data2!
    
    Now at batch 0<->32...
    42 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 90.64999999999999%! 
    
    There are 26 unmatched samples and 252 matched samples in data1!
    
    53 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.41000000000001%! 
    
    There are 17 unmatched samples and 287 matched samples in data2!
    
    Now at batch 0<->33...
    49 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 92.81%! 
    
    There are 20 unmatched samples and 258 matched samples in data1!
    
    56 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.08%! 
    
    There are 18 unmatched samples and 286 matched samples in data2!
    
    Now at batch 0<->34...
    45 cells from Modality X are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 94.24%! 
    
    There are 16 unmatched samples and 262 matched samples in data1!
    
    53 cells from Modality Y are unmatched in Phase I and are realigned in Phase II.
    iterative unbalanced optimal transport converges after 2 iterations with cell-type matching accuracy 95.75%! 
    
    There are 13 unmatched samples and 293 matched samples in data2!
    


.. code:: ipython3

    rna_source_ct_array = combined_adata.obs['cell_type']
    protein_source_ct_array = protein_adata_filter.obs['cell_type']
    rna_aligned_protein, rna_predict_ct_array, protein_aligned_rna, protein_predict_ct_array = cellink.synchronize_imputed_to_initial()

.. code:: ipython3

    protein_aligned_rna = protein_aligned_rna[protein_source_ct_array == protein_predict_ct_array, :]
    ct_protein_aligned = protein_source_ct_array[protein_source_ct_array == protein_predict_ct_array]
    rna_aligned_protein = rna_aligned_protein[rna_source_ct_array == rna_predict_ct_array, :]
    ct_rna_aligned = rna_source_ct_array[rna_source_ct_array == rna_predict_ct_array]

.. code:: ipython3

    rna_matched_cellids = rna_source_ct_array == rna_predict_ct_array
    protein_matched_cellids = protein_source_ct_array == protein_predict_ct_array

.. code:: ipython3

    cellink.visualize_integration(ann1_full_batch=combined_adata, ann2_full_batch=protein_adata_filter, arr2_imputed=protein_aligned_rna,
                                  datatype= ['scRNA', 'CODEX'], matched_cellids=protein_matched_cellids, direction=1)



.. image:: HPAP023_analysis_files/HPAP023_analysis_17_0.png


.. code:: ipython3

    cellink.visualize_integration(ann1_full_batch=protein_adata_filter, ann2_full_batch=combined_adata, arr2_imputed=rna_aligned_protein,
                                  datatype= ['CODEX', 'scRNA'], matched_cellids=rna_matched_cellids, direction=2)



.. image:: HPAP023_analysis_files/HPAP023_analysis_18_0.png

