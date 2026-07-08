# ==============================================================================
# scripts/10_graph_embedding_validation.py
# ------------------------------------------------------------------------------
# [EN] Week 7: Graph Embedding Validation using STRING and CTD databases.
# [ZH] 第 7 周：利用 STRING 和 CTD 数据库进行图嵌入验证。
# ==============================================================================

import os
import gzip
import shutil
import urllib.request
import pandas as pd
import numpy as np
import networkx as nx
from node2vec import Node2Vec
import umap
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity

# Setup directories
os.makedirs("data_raw/string", exist_ok=True)
os.makedirs("data_raw/ctd", exist_ok=True)
os.makedirs("data_clean", exist_ok=True)
os.makedirs("Week7", exist_ok=True)
os.makedirs("results", exist_ok=True)

# ==============================================================================
# SECTION 1: Download Datasets | 第 1 部分：下载数据库文件
# ==============================================================================

string_links_url = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
string_alias_url = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
ctd_gda_url = "https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz"

string_links_dest = "data_raw/string/9606.protein.links.v12.0.txt.gz"
string_alias_dest = "data_raw/string/9606.protein.aliases.v12.0.txt.gz"
ctd_gda_dest = "data_raw/ctd/CTD_genes_diseases.tsv.gz"

def download_file(url, dest):
    if os.path.exists(dest):
        print(f"File already exists: {dest}")
        return
    print(f"Downloading {url} to {dest}...")
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36'
    }
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req) as response, open(dest, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print("Download completed successfully.")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        raise e

download_file(string_links_url, string_links_dest)
download_file(string_alias_url, string_alias_dest)
download_file(ctd_gda_url, ctd_gda_dest)

# ==============================================================================
# SECTION 2: Load and Map STRING Gene Aliases (Priority-based) | 第 2 部分：映射基因别名
# ==============================================================================

def load_string_aliases(alias_path):
    print("Parsing STRING protein-to-gene mapping...")
    protein_to_gene = {}
    with gzip.open(alias_path, 'rt', encoding='utf-8') as f:
        # Read header
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                prot_id, alias, source = parts[0], parts[1], parts[2]
                
                # Check for high-priority sources
                if source in ['Ensembl_HGNC_symbol', 'BioMart_HUGO', 'UniProt_GN_Name', 'Ensembl_HGNC']:
                    # High priority always overwrites or initializes (Score 2)
                    protein_to_gene[prot_id] = (alias, 2)
                elif source in ['Ensembl_EntrezGene', 'Ensembl_WikiGene', 'KEGG_NAME']:
                    # Medium priority only writes if not already mapped or if mapped with low priority (Score 1)
                    if prot_id not in protein_to_gene or protein_to_gene[prot_id][1] < 2:
                        protein_to_gene[prot_id] = (alias, 1)

    # Convert to simple mapping dictionary
    mapping = {k: v[0] for k, v in protein_to_gene.items()}
    print(f"Mapped {len(mapping)} protein IDs to Gene Symbols.")
    return mapping

protein_to_gene = load_string_aliases(string_alias_dest)

# ==============================================================================
# SECTION 3: Build Heterogeneous Graph G | 第 3 部分：构建全基因组异构网络
# ==============================================================================

G = nx.Graph()

# 3.1 Add STRING PPI Edges (High Confidence)
print("Loading STRING PPI links...")
high_conf_edges = 0
with gzip.open(string_links_dest, 'rt', encoding='utf-8') as f:
    # Read header
    next(f)
    for line in f:
        parts = line.strip().split(' ')
        if len(parts) >= 3:
            p1, p2, score = parts[0], parts[1], int(parts[2])
            if score >= 700:  # High confidence interaction
                g1 = protein_to_gene.get(p1)
                g2 = protein_to_gene.get(p2)
                if g1 and g2:
                    G.add_edge(g1, g2, weight=score / 1000.0, edge_type='PPI')
                    high_conf_edges += 1

print(f"Added {high_conf_edges} high-confidence PPI edges to the graph.")

# 3.2 Add CTD curated Gene-Disease Associations (GDAs)
print("Loading CTD gene-disease associations...")
gda_df = pd.read_csv(ctd_gda_dest, sep='\t', comment='#', header=None,
                     names=['GeneSymbol', 'GeneID', 'DiseaseName', 'DiseaseID', 
                            'DirectEvidence', 'InferenceChemicalName', 
                            'InferenceScore', 'OmimIDs', 'PubMedIDs'])

# Filter for three target diseases of interest
target_diseases = ['Leukemia, Myeloid, Acute', 'Lymphoma, Large B-Cell, Diffuse', 'Multiple Myeloma']
gda_df = gda_df[gda_df['DiseaseName'].isin(target_diseases)]

# Keep GDA rows with direct curated evidence or high inference score >= 30.0
gda_df = gda_df[(gda_df['InferenceScore'] >= 30.0) | gda_df['DirectEvidence'].notna()]

added_gda_edges = 0
for idx, row in gda_df.iterrows():
    gene = str(row['GeneSymbol']).strip()
    disease_id = str(row['DiseaseID']).strip()
    disease_name = str(row['DiseaseName']).strip()
    
    # We add disease nodes prefixed with 'Disease:' to separate namespaces
    disease_node = f"Disease:{disease_id}"
    
    # Add gene to graph if not present, and link to disease
    if gene and disease_id and gene in G:
        G.add_node(disease_node, label='Disease', name=disease_name)
        ref_score = row['InferenceScore']
        weight = float(ref_score) if pd.notna(ref_score) else 1.0
        G.add_edge(gene, disease_node, weight=weight, edge_type='GDA')
        added_gda_edges += 1

print(f"Added {added_gda_edges} curated Gene-Disease Association edges.")
print(f"Total graph size: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

# ==============================================================================
# SECTION 4: Train Node2Vec Embedding Model | 第 4 部分：训练 Node2Vec 节点嵌入
# ==============================================================================

print("Generating random walks and training Node2Vec (this may take 1-2 minutes)...")
# walk_length=30, num_walks=5 for optimized execution in this environment
node2vec_model = Node2Vec(
    G, 
    dimensions=128, 
    walk_length=30, 
    num_walks=5, 
    p=1.0, 
    q=0.5, 
    workers=4
)

print("Fitting Word2Vec model...")
model = node2vec_model.fit(window=10, min_count=1, batch_words=4)

# Save embeddings in text format (immune to decoding issues)
model.wv.save_word2vec_format("data_clean/node2vec_embeddings.txt", binary=False)
print("Embeddings saved to data_clean/node2vec_embeddings.txt")

# ==============================================================================
# SECTION 5: Proximity Analysis & Re-ranking | 第 5 部分：图嵌入空间邻近度计算与重排序
# ==============================================================================

# Define Candidates (Top 10 from Week 5/6)
candidates = ["IL1R2", "VNN2", "SLC15A3", "EGFL7", "CXCR1", "CXCR2", "CD14", "GZMB", "FPR2", "CMTM2"]

# Define Known AML Targets
known_targets = ["CD33", "IL3RA", "FLT3", "BCL2", "IDH1", "IDH2", "KIT", "CD47", "HAVCR2", "WT1", "KMT2A", "TP53", "CD38", "NPM1"]

# Filter for those present in the trained model vocabulary
valid_candidates = [c for c in candidates if c in model.wv]
valid_known = [k for k in known_targets if k in model.wv]

print(f"\nValid Candidates in model vocabulary: {valid_candidates}")
print(f"Valid Known Targets in model vocabulary: {valid_known}")

# Compute centroid of known targets
known_vectors = np.array([model.wv[k] for k in valid_known])
centroid = np.mean(known_vectors, axis=0).reshape(1, -1)

# Calculate cosine similarities
proximity_results = []
for gene in valid_candidates:
    gene_vec = model.wv[gene].reshape(1, -1)
    sim = cosine_similarity(gene_vec, centroid)[0][0]
    
    proximity_results.append({
        'gene_symbol': gene,
        'graph_cosine_similarity': sim
    })

proximity_df = pd.DataFrame(proximity_results).sort_values(by='graph_cosine_similarity', ascending=False)
proximity_df['graph_rank'] = range(1, len(proximity_df) + 1)

# Save results
proximity_df.to_csv("data_clean/target_graph_proximity_rankings.csv", index=False)
proximity_df.to_csv("Week7/target_graph_proximity_rankings.csv", index=False)

print("\n--- Graph Embedding Proximity Rankings ---")
print(proximity_df)

# ==============================================================================
# SECTION 6: UMAP Dimensionality Reduction & Plotting | 第 6 部分：UMAP 降维与散点可视化
# ==============================================================================

print("\nRunning UMAP dimensionality reduction on all gene embeddings...")
# Get all gene nodes in the model vocabulary (excluding disease nodes prefixed with 'Disease:')
all_genes_in_vocab = [node for node in model.wv.index_to_key if not node.startswith("Disease:")]
all_gene_vectors = np.array([model.wv[node] for node in all_genes_in_vocab])

reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_2d = reducer.fit_transform(all_gene_vectors)

# Create plotting dataframe
plot_df = pd.DataFrame({
    'gene_symbol': all_genes_in_vocab,
    'UMAP1': embedding_2d[:, 0],
    'UMAP2': embedding_2d[:, 1],
    'Category': 'Background Genes'
})

# Annotate Categories
plot_df.loc[plot_df['gene_symbol'].isin(valid_known), 'Category'] = 'Known AML Targets'
plot_df.loc[plot_df['gene_symbol'].isin(valid_candidates), 'Category'] = 'Our Candidates'

# Plot
plt.figure(figsize=(10, 8))

# Background
bg = plot_df[plot_df['Category'] == 'Background Genes']
plt.scatter(bg['UMAP1'], bg['UMAP2'], c='#D9D9D9', alpha=0.4, s=6, label='Background Genes')

# Known Targets
known = plot_df[plot_df['Category'] == 'Known AML Targets']
plt.scatter(known['UMAP1'], known['UMAP2'], c='#E41A1C', alpha=0.9, s=80, marker='*', label='Known AML Targets', edgecolors='black')

# Our Candidates
cands = plot_df[plot_df['Category'] == 'Our Candidates']
plt.scatter(cands['UMAP1'], cands['UMAP2'], c='#377EB8', alpha=0.9, s=70, marker='D', label='Our Candidates', edgecolors='black')

# Label Candidates
for idx, row in cands.iterrows():
    plt.annotate(
        row['gene_symbol'], 
        (row['UMAP1'], row['UMAP2']),
        textcoords="offset points", 
        xytext=(0, 6), 
        ha='center', 
        fontsize=9, 
        fontweight='bold',
        bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.6, ec="orange")
    )

plt.title("UMAP Proximity of Target Candidates in Human Interactome Network", fontsize=12, fontweight='bold')
plt.xlabel("UMAP Dimension 1", fontsize=10)
plt.ylabel("UMAP Dimension 2", fontsize=10)
plt.legend(loc='best', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig("Week7/umap_embedding_proximity_plot.png", dpi=300)
plt.savefig("results/umap_embedding_proximity_plot.png", dpi=300)
plt.close()

print("UMAP proximity plots saved successfully to Week7/ and results/")
print("Graph validation pipeline completed successfully!")
