library(tidyverse)
library(Matrix)
library(Seurat)
sanes_count_matrix <- data.table::fread('GSE149715_MouseAC_count_matrix.csv.gz') 
rownames(sanes_count_matrix) <- sanes_count_matrix$V1
sanes_count_matrix <- sanes_count_matrix[,-1]
sanes_count_matrix <-Matrix(as.matrix(sanes_count_matrix, sparse = T) )
save(sanes_count_matrix, file = 'references/sanes_amacrine_sparse_matrix.Rdata')
sanes_meta_data <- read_tsv('sanes_amacrine_metadata.tsv.gz', col_names = c('sample_accession', 'BioSample', 'mouse', 'path'))
sanes_metadata_full <- tibble(Barcode = colnames(sanes_count_matrix)) %>% mutate(mouse = str_extract(Barcode, 'MouseACS\\d+'),
                                                                                 mouse = as.factor(mouse)) %>% 
  inner_join(sanes_meta_data) %>% as.data.frame
rownames(sanes_metadata_full) <- sanes_metadata_full$Barcode

seu <- CreateSeuratObject(sanes_count_matrix, meta.data = sanes_metadata_full )
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
seu <- ScaleData(seu, vars.to.regress = 'mouse')
# #save(seu, file = 'testing/sanes_amacrine_seu_pp.Rdata')
# load('/data/swamyvs/scEiad_subcelltype/testing/sanes_amacrine_seu_pp.Rdata')
seu <- RunPCA(seu, npcs = 40,features = VariableFeatures(object = seu))
seu_pca_data = Embeddings(seu, 'pca') %>% as.data.frame %>% rownames_to_column('Barcode')
write_csv(seu_pca_data, 'data/sanes_amacrine_preproccessed.csv.gz')
sanes_amacrine_labels <- read_csv('/data/swamyvs/scEiad_subcelltype/references/sanes_amacrine_subcluster_labels.csv', 
                                  skip=2, col_names = c('Barcode', 'Group'))
seu@meta.data %>% select(Barcode) %>% inner_join(sanes_amacrine_labels) %>% dplyr::rename(sanes_bc=Barcode) %>% 
  mutate(Barcode = 0:(nrow(.)-1)) %>% write_csv('data/true_labels/sanes_amacrine_labels.csv.gz')


pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 40,features = VariableFeatures(object = pbmc))
pca_data =  Embeddings(pbmc, 'pca') %>% as.data.frame %>% rownames_to_column('Barcode')
write_csv(pca_data, 'data/pbmc_preproccessed.csv.gz')
pbmc_labels <- read_csv('data/pbmc_celltypes.csv.gz')
pbmc@meta.data %>% rownames_to_column('Barcode') %>%  select(Barcode) %>% left_join(pbmc_labels) %>% 
  dplyr::rename(sanes_bc=Barcode, Group = cellType) %>% 
  mutate(Barcode = 0:(nrow(.)-1), Group = replace_na(Group, 'MSG')) %>% 
  write_csv('data/true_labels/pbmc_labels.csv.gz')


