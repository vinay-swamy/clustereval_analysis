import numpy as np
import re 

rule all:
    input:
        clustering = expand('cluster_out/ds-{countMat}_alg-{alg}_knn-{k}_.csv.gz', k=list(range(5,101,5)),
                                countMat = ['sanes_amacrine','pbmc', 'pbmc_sim_20groups-V1', 'pbmc_sim_10groups-V1'], 
                                alg = ['louvain', 'leiden', 'louvainPrune', 'leidenPrune']),
        metrics = expand('cluster_metrics/ds-{countMat}_alg-{alg}_knn-{k}_.csv.gz', k=list(range(5,101,5)),
                                countMat = ['sanes_amacrine','pbmc', 'pbmc_sim_20groups-V1', 'pbmc_sim_10groups-V1'], 
                                alg = ['louvain', 'leiden','louvainPrune', 'leidenPrune']),
        pairwise = expand('cluster_metrics/ds_pairwise-{countMat}-exp_alg-{alg}_metrics.pickle',
                        countMat = ['sanes_amacrine','pbmc', 'pbmc_sim_20groups-V1', 'pbmc_sim_10groups-V1'],
                        alg = ['louvain', 'leiden','louvainPrune', 'leidenPrune'])

rule make_count_matrices:
    output:
        'data/count_mats/sanes_amacrine_preproccessed.csv.gz',
        'data/count_mats/pbmc_preproccessed.csv.gz',
    shell:
        '''

        wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
        wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149715/suppl/GSE149715%5FMouseAC%5Fcount%5Fmatrix%2Ecsv%2Egz
        tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
        module load R/4.0.3
        Rscript scripts/make_raw_data.R
        rm pbmc3k_filtered_gene_bc_matrices.tar.gz
        rm filtered_gene_bc_matrices/ -rf
        rm GSE149715_MouseAC_count_matrix.csv.gz

        '''


rule run_parameter_exp_clustering:
    input: 
        mat = 'data/count_mats/{countMat}_preproccessed.csv.gz'
    output:
        clustering = 'cluster_out/ds-{countMat}_alg-{alg}_knn-{k}_.csv.gz',
        metrics = 'cluster_metrics/ds-{countMat}_alg-{alg}_knn-{k}_.csv.gz'

    run:
        import clustereval as ce
        import pandas as pd 
        
        if wildcards.alg in ['louvainPrune', 'leidenPrune'] :
            
            metrics,labels, perturbations = ce.cluster.run_full_experiment(reduction = pd.read_csv(input.mat), 
                                            alg = re.sub('Prune', '', wildcards.alg), 
                                            k=int(wildcards.k),
                                            n_perturbations=100,
                                            edge_permut_frac=.05,
                                            local_pruning_dist_threshold=3,
                                            global_pruning_jac_threshold = 'median',
                                            min_cluster_size=25,
                                            experiment_name=f'experiment2-{wildcards.countMat}-{str(wildcards.k)}-{wildcards.alg}'
                                            )
        else:
            metrics,labels, perturbations = ce.cluster.run_full_experiment(reduction = pd.read_csv(input.mat), 
                                alg = wildcards.alg, 
                                k=int(wildcards.k),
                                n_perturbations=100,
                                edge_permut_frac=.05,
                                min_cluster_size=25,
                                experiment_name=f'experiment2-{wildcards.countMat}-{str(wildcards.k)}-{wildcards.alg}'
                                )

        
        labels.to_csv(output.clustering,index=False,header=False)
        metrics.to_csv(output.metrics)



rule calculate_pairwise_metrics:
    input: 
        expand('cluster_out/ds-{{countMat}}_alg-{{alg}}_knn-{k}_.csv.gz', 
                k=list(range(5,101,5)) )
    output:
        results = 'cluster_metrics/ds_pairwise-{countMat}-exp_alg-{alg}_metrics.pickle'
    params:
        glob_string= lambda wildcards: f'cluster_out/ds-{wildcards.countMat}_alg-{wildcards.alg}_knn-*_.csv.gz'
    run:
        import clustereval as ce
        import numpy as np 
        import pandas as pd
        import pickle
        exp_results = ce.metrics.pairwise_metric_calculation(params.glob_string, 96)
        with open(output.results, 'wb+') as outfile:
            exp_results = [{'exp_param': i.exp_param, 
                            'cluster_ids':i.cluster_ids, 
                            'stability_scores': i.stability_scores, 
                            'purity_scores': i.purity_scores} for i in exp_results  ]
            pickle.dump(exp_results, outfile)
        
