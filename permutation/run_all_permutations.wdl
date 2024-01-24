# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Scatter germline:somatic permutation testing over one or more sets of empirical data


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/germline-somatic-exploration-2023/main/permutation/scatter_permutation_single_strata.wdl" as Permute


workflow ScatterPermutations {
  input {
    # Empirically observed overlaps
    Array[File] observed_overlaps_tsvs
    Array[File] strata_gene_count_tsvs
    Array[File] expression_quantile_gene_counts_tsvs
    Array[String] overlap_prefixes

    # Data required for each permutation
    File uniform_weights
    File coding_nonsyn_weights
    File coding_gwas_weights
    File noncoding_gwas_weights
    File somatic_noncoding_weights
    File eligible_gene_symbols
    File expression_quantiles

    # Reference files for assessing overlap
    File cellchat_tsv
    File ppi_tsv
    File complexes_tsv

    # Certain scripts supplied as inputs to avoid needing to build a special
    # Docker image specifically for this analysis
    File shuffle_script
    File find_pairs_script
    File analysis_script

    # Parallelization control
    Int n_shards = 100
    Int n_perm_per_shard = 100
    Int n_cpu_per_shard = 2
    Float? mem_gb_per_shard
    Int disk_gb_per_shard = 30

    # Analysis control
    Int analysis_n_cpu = 4
    Float? analysis_mem_gb
    Int analysis_disk_gb = 30
    
    String docker = "vanallenlab/g2c_pipeline:r_argparse"
  }

  Int n_sets = length(observed_overlaps_tsvs)

  scatter ( i in range(n_sets) ) {
    call Permute.RunPermutations as RunPermutations {
      input:
      observed_overlaps_tsv = observed_overlaps_tsvs[i],
      output_prefix = overlap_prefixes[i],
      strata_gene_counts_tsv = strata_gene_count_tsvs[i],
      uniform_weights = uniform_weights,
      coding_nonsyn_weights = coding_nonsyn_weights,
      coding_gwas_weights = coding_gwas_weights,
      noncoding_gwas_weights = noncoding_gwas_weights,
      somatic_noncoding_weights = somatic_noncoding_weights,
      eligible_gene_symbols = eligible_gene_symbols,
      expression_quantiles = expression_quantiles,
      expression_quantile_gene_counts_tsv = expression_quantile_gene_counts_tsvs[i],
      cellchat_tsv = cellchat_tsv,
      ppi_tsv = ppi_tsv,
      complexes_tsv = complexes_tsv,
      shuffle_script = shuffle_script,
      find_pairs_script = find_pairs_script,
      analysis_script = analysis_script,
      n_shards = n_shards,
      n_perm_per_shard = n_perm_per_shard,
      n_cpu_per_shard = n_cpu_per_shard,
      mem_gb_per_shard = mem_gb_per_shard,
      disk_gb_per_shard = disk_gb_per_shard,
      analysis_n_cpu = analysis_n_cpu,
      analysis_mem_gb = analysis_mem_gb,
      analysis_disk_gb = analysis_disk_gb,
      docker = docker
    }
  }

  call Permute.BundleOutputs as BundleOutputs {
    input:
      tarballs = RunPermutations.results,
      output_prefix = "all_permutation_results",
      docker = docker
  }

  output {
    File results_tarball = BundleOutputs.results_tarball
  }
}

