import clustereval
from clustereval import cluster 
import louvain 
import pandas as pd 
def louvain_clustering(reduction, res, k, perturb = False):
    
    clu_obj = cluster.Cluster(data=reduction, knn=k,  nthreads=1)
    clu_obj.buildNeighborGraph(nn_space='l2', ef_construction=150,
                            local_pruning=False, global_pruning=False, jac_std_global='median')
    if perturb:
        clu_obj.run_perturbation()
    labels = clu_obj.run_louvain(
        vertex_partition_method=louvain.RBConfigurationVertexPartition,
        resolution=res,
        jac_weighted_edges='weight'
    )
    labels_corrected = clu_obj.merge_singletons(labels, 10)
    outdf = pd.DataFrame(
        {"Barcode": list(reduction.index), 'labels': labels_corrected})

    return outdf 

reduction = pd.read_csv('/data/swamyvs/scEiad_subcelltype/testing/amacrin_mm_scvi_dim.csv', index_col=0)
print('BEGIN CLUSTER ref')
ref_df = louvain_clustering(reduction, 1.0, 30, perturb=False).sort_values('labels').assign(
    labels=lambda x: 'clu_' + x.labels.astype(str)).to_dict('list')
print('BEGIN CLUSTER other ')
exp_dfs = [louvain_clustering(reduction, 1.0, i, perturb=False).sort_values('labels').assign(
    labels=lambda x: 'clu_' + x.labels.astype(str)).to_dict(
    'list') for i in range(30, 60, 5)]
cluster_names = [f'louvain-1.0-{str(i)}' for i in range(30, 60, 5)]
print('BEGIN ONEWAY')
clustereval.calc_metrics.oneway_metric_calculation(ref_df, exp_dfs, 'test', 'test_out.csv')

print('BEGIN PAIRWISE')
res = clustereval.calc_metrics.pairwise_metric_calculation_frommem(
    exp_dfs,cluster_names, 4)
print(res)
