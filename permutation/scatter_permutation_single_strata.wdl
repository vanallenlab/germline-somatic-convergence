# Germline:Somatic Exploratory Pilot Analysis 2023
# Copyright (c) 2023 Sam Hoffman, Ryan Collins, and the Van Allen Lab @ Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL 2.0 (see LICENSE)

# Parallelizable permutation of gene lists to assess germline:somatic overlap


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utils


workflow RunPermutations {
  input {
    # Empirically observed overlaps
    File observed_overlaps_tsv
    File expression_quantile_gene_counts_tsv
    String output_prefix

    # Data required for each permutation
    File strata_gene_counts_tsv
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

  Int total_n_perms = n_shards * n_perm_per_shard

  scatter ( shard_number in range(n_shards) ) {

    # One set of permutations using uniform weights for all strata
    call PermuteOverlaps as PermuteUniform {
      input:
        strata_gene_counts_tsv = strata_gene_counts_tsv,
        coding_nonsyn_weights = uniform_weights,
        coding_gwas_weights = uniform_weights,
        noncoding_gwas_weights = uniform_weights,
        somatic_noncoding_weights = uniform_weights,
        eligible_gene_symbols = eligible_gene_symbols,
        cellchat_tsv = cellchat_tsv,
        ppi_tsv = ppi_tsv,
        complexes_tsv = complexes_tsv,
        shuffle_script = shuffle_script,
        find_pairs_script = find_pairs_script,
        shard_number = shard_number,
        n_perm_per_shard = n_perm_per_shard,
        docker = docker,
        n_cpu = n_cpu_per_shard,
        mem_gb = mem_gb_per_shard,
        disk_gb = disk_gb_per_shard
    }

    # One set of permutations using custom prior weights
    call PermuteOverlaps as PermuteBayesian {
      input:
        strata_gene_counts_tsv = strata_gene_counts_tsv,
        coding_nonsyn_weights = coding_nonsyn_weights,
        coding_gwas_weights = coding_gwas_weights,
        noncoding_gwas_weights = noncoding_gwas_weights,
        somatic_noncoding_weights = somatic_noncoding_weights,
        eligible_gene_symbols = eligible_gene_symbols,
        cellchat_tsv = cellchat_tsv,
        ppi_tsv = ppi_tsv,
        complexes_tsv = complexes_tsv,
        shuffle_script = shuffle_script,
        find_pairs_script = find_pairs_script,
        shard_number = shard_number,
        n_perm_per_shard = n_perm_per_shard,
        docker = docker,
        n_cpu = n_cpu_per_shard,
        mem_gb = mem_gb_per_shard,
        disk_gb = disk_gb_per_shard
    }

    # One set of permutations with tissue-specific expression quantile matching
    call PermuteOverlaps as PermuteExpression {
      input:
        strata_gene_counts_tsv = strata_gene_counts_tsv,
        coding_nonsyn_weights = uniform_weights,
        coding_gwas_weights = uniform_weights,
        noncoding_gwas_weights = uniform_weights,
        somatic_noncoding_weights = uniform_weights,
        eligible_gene_symbols = eligible_gene_symbols,
        expression_quantiles = expression_quantiles,
        expression_quantile_gene_counts_tsv = expression_quantile_gene_counts_tsv,
        cellchat_tsv = cellchat_tsv,
        ppi_tsv = ppi_tsv,
        complexes_tsv = complexes_tsv,
        shuffle_script = shuffle_script,
        find_pairs_script = find_pairs_script,
        shard_number = shard_number,
        n_perm_per_shard = n_perm_per_shard,
        docker = docker,
        n_cpu = n_cpu_per_shard,
        mem_gb = mem_gb_per_shard,
        disk_gb = disk_gb_per_shard
    }

    # One set of permutations combining expression quantile-matching and custom prior weights
    call PermuteOverlaps as PermuteComposite {
      input:
        strata_gene_counts_tsv = strata_gene_counts_tsv,
        coding_nonsyn_weights = coding_nonsyn_weights,
        coding_gwas_weights = coding_gwas_weights,
        noncoding_gwas_weights = noncoding_gwas_weights,
        somatic_noncoding_weights = somatic_noncoding_weights,
        eligible_gene_symbols = eligible_gene_symbols,
        expression_quantiles = expression_quantiles,
        expression_quantile_gene_counts_tsv = expression_quantile_gene_counts_tsv,
        cellchat_tsv = cellchat_tsv,
        ppi_tsv = ppi_tsv,
        complexes_tsv = complexes_tsv,
        shuffle_script = shuffle_script,
        find_pairs_script = find_pairs_script,
        shard_number = shard_number,
        n_perm_per_shard = n_perm_per_shard,
        docker = docker,
        n_cpu = n_cpu_per_shard,
        mem_gb = mem_gb_per_shard,
        disk_gb = disk_gb_per_shard
    }
  }

  # Concatenate results from uniformly weighted permutations
  call Utils.ConcatTextFiles as ConcatUniform {
    input:
      shards = PermuteUniform.results_tsv,
      compression_command = "gzip -c",
      output_filename = "germ_som_convergence.permutation_results.uniform_weighting.n" + total_n_perms + ".txt.gz",
      docker = docker
  }

  # Compare uniformly permuted results to empirically observed results
  call ComparePermutedAndEmpirical as AnalyzeUniform {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatUniform.merged_file,
      cancers = PermuteUniform.cancers[0],
      analysis_script = analysis_script,
      output_prefix = "uniform_permutation",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from Bayesian weighted permutations
  call Utils.ConcatTextFiles as ConcatBayesian {
    input:
      shards = PermuteBayesian.results_tsv,
      compression_command = "gzip -c",
      output_filename = "germ_som_convergence.permutation_results.nonuniform_weighting.n" + total_n_perms + ".txt.gz",
      docker = docker
  }

  # Compare Bayesian permuted results to empirically observed results
  call ComparePermutedAndEmpirical as AnalyzeBayesian {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatBayesian.merged_file,
      cancers = PermuteBayesian.cancers[0],
      analysis_script = analysis_script,
      output_prefix = "nonuniform_permutation",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from expression quantile-matched permutations
  call Utils.ConcatTextFiles as ConcatExpression {
    input:
      shards = PermuteExpression.results_tsv,
      compression_command = "gzip -c",
      output_filename = "germ_som_convergence.permutation_results.expression_weighting.n" + total_n_perms + ".txt.gz",
      docker = docker
  }

  # Compare expression quantile-matched permuted results to empirically observed results
  call ComparePermutedAndEmpirical as AnalyzeExpression {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatExpression.merged_file,
      cancers = PermuteExpression.cancers[0],
      analysis_script = analysis_script,
      output_prefix = "expression_permutation",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from composite weighted permutations
  call Utils.ConcatTextFiles as ConcatComposite {
    input:
      shards = PermuteComposite.results_tsv,
      compression_command = "gzip -c",
      output_filename = "germ_som_convergence.permutation_results.composite_weighting.n" + total_n_perms + ".txt.gz",
      docker = docker
  }

  # Compare composite permuted results to empirically observed results
  call ComparePermutedAndEmpirical as AnalyzeComposite {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatComposite.merged_file,
      cancers = PermuteComposite.cancers[0],
      analysis_script = analysis_script,
      output_prefix = "composite_permutation",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Bundle all results as a single tarball for returning to parent workflow
  call BundleOutputs {
    input:
      tarballs = [AnalyzeUniform.results_tarball, AnalyzeBayesian.results_tarball, AnalyzeExpression.results_tarball, AnalyzeComposite.results_tarball],
      output_prefix = output_prefix,
      docker = docker
  }

  output {
    File results = BundleOutputs.results_tarball
  }
}


task PermuteOverlaps {
  input {
    File strata_gene_counts_tsv
    File coding_nonsyn_weights
    File coding_gwas_weights
    File noncoding_gwas_weights
    File somatic_noncoding_weights
    File eligible_gene_symbols
    File? expression_quantiles
    File? expression_quantile_gene_counts_tsv

    File cellchat_tsv
    File ppi_tsv
    File complexes_tsv

    File shuffle_script
    File find_pairs_script

    Int shard_number
    Int n_perm_per_shard

    String docker
    Int n_cpu = 2
    Float? mem_gb
    Int disk_gb = 30
  }

  String outfile_name = "perm_results." + shard_number + ".tsv"

  command <<<
    set -eu -o pipefail

    # Get list of cancer types
    cut -f1 ~{strata_gene_counts_tsv} | sort -V | uniq > cancers.list

    # Run multiple permutations in each task
    for i in $( seq 1 ~{n_perm_per_shard} ); do

      echo -e "Starting permutation $i of ~{n_perm_per_shard}..."

      # Reset directory tree for each permutation
      if [ -e gene_lists ]; then
        rm -rf gene_lists
      fi
      mkdir gene_lists
      for origin in germline somatic; do
        for context in coding noncoding; do
          mkdir gene_lists/${origin}_${context}
        done
      done
      mkdir \
        gene_lists/germline_coding_cosmic \
        gene_lists/germline_coding_gwas
      if [ -e overlaps ]; then
        rm -rf overlaps
      fi
      mkdir overlaps
      if [ -e permutation_results.tmp.tsv ]; then
        rm permutation_results.tmp.tsv
      fi

      # Permute each strata with related but distinct seeds
      while read cancer origin context n_genes; do

        # Define strata-specific weights
        case "${origin}_${context}" in
          germline_coding_gwas)
            weights=~{coding_gwas_weights}
            # For coding GWAS hits, by design we want to exclude genes that were
            # already sampled from the coding germline COSMIC set
            fgrep \
              -wvf gene_lists/germline_coding_cosmic/$cancer.germline.coding_cosmic.genes.list \
              ~{eligible_gene_symbols} \
            > tmp.coding_gwas.elig_genes.list
            elig_list=tmp.coding_gwas.elig_genes.list
            ;;
          germline_noncoding)
            weights=~{noncoding_gwas_weights}
            elig_list=~{eligible_gene_symbols}
            ;;
          somatic_noncoding)
            weights=~{somatic_noncoding_weights}
            elig_list=~{eligible_gene_symbols}
            ;;
          *)
            weights=~{coding_nonsyn_weights}
            elig_list=~{eligible_gene_symbols}
            ;;
        esac

        # Set seed string
        # Note: we take the first and last characters of each cancer, origin, 
        # and context because of total seed string length limits in our R script
        sp1="~{shard_number}"
        sp2="$( echo $cancer | cut -c1 )${cancer: -1}"
        sp3="$( echo $origin | cut -c1 )${origin: -1}"
        sp4="$( echo $context | cut -c1 )${context: -1}"
        sp5="$i"
        seed="${sp1}${sp2}${sp3}${sp4}${sp5}"

        if [ ~{defined(expression_quantiles)} == "true" ] && \
           [ ~{defined(expression_quantile_gene_counts_tsv)} == "true"]; then
          # Shuffle & sample each expression quantile separately, if optioned
          while read q n; do
            if [ $n -eq 0 ]; then
              continue
            fi
            # Define list of genes in this expression quantile in this tissue
            awk -v FS="\t" -v cancer=$cancer -v q=$q \
              '{ if ($2==q && $3==cancer) print $1 }' \
              ~{expression_quantiles} \
            | fgrep -xf $elig_list \
            | sort -V > $cancer.gex_q$q.genes.list
            # Shuffle genes, restrict to eligible gene symbols, and sample desired number
            Rscript ~{shuffle_script} \
              --tsv-in $weights \
              --seed "$seed" \
              --eligible-genes $cancer.gex_q$q.genes.list \
              --outfile tmp.shuffled.genes.list
            # Note: need to break this up to avoid Rcript throwing SIGPIPE error
            head -n $n tmp.shuffled.genes.list \
            >> gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          done < <( awk -v FS="\t" -v OFS="\t" -v cancer=$cancer -v origin=$origin -v context=$context \
                      '{ if ($1==cancer && $3==origin && $4==context) print $2, $5 }' \
                      ~{expression_quantile_gene_counts_tsv} )
        else
          if [ $n_genes -gt 0 ]; then
            # Shuffle genes, restrict to eligible gene symbols, and sample desired number
            Rscript ~{shuffle_script} \
              --tsv-in $weights \
              --seed "$seed" \
              --eligible-genes $elig_list \
              --outfile tmp.shuffled.genes.list
            # Note: need to break this up to avoid Rcript throwing SIGPIPE error
            head -n $n_genes tmp.shuffled.genes.list \
            > gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          else
            touch gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          fi
        fi

      done < ~{strata_gene_counts_tsv}

      # Find pairs per cancer type
      while read cancer; do

        # Merge germline coding COSMIC & coding GWAS genes
        cat \
          gene_lists/germline_coding_cosmic/$cancer.germline.coding_cosmic.genes.list \
          gene_lists/germline_coding_gwas/$cancer.germline.coding_gwas.genes.list \
        | sort -V | uniq \
        > gene_lists/germline_coding/$cancer.germline.coding.genes.list

        # Find all pairs & summarize data as simple .tsv line
        for germ_context in coding noncoding; do
          for som_context in coding noncoding; do
            python ~{find_pairs_script} \
              --germline gene_lists/germline_$germ_context/$cancer.germline.$germ_context.genes.list \
              --somatic gene_lists/somatic_$som_context/$cancer.somatic.$som_context.genes.list \
              --cellchat-db ~{cellchat_tsv} \
              --ppi-db ~{ppi_tsv} \
              --protein-complexes ~{complexes_tsv} \
              --report-counts
          done
        done \
        | paste -s \
        | awk -v OFS="\t" -v cancer=$cancer -v major_idx=~{shard_number} -v minor_idx=$i \
          '{ print cancer, major_idx, minor_idx, 
                   $1+$6+$11+$16, $2+$7+$12+$17, $3+$8+$13+$18, \
                   $4+$9+$14+$19, $5+$10+$15+$20, $0 }' \
        >> permutation_results.tmp.tsv

      done < cancers.list

      # Add row for pan-cancer analysis (column-wise sum of all cancers)
      cut -f4- permutation_results.tmp.tsv \
      | awk -v OFS="\t" -v cancer="all" -v major_idx=~{shard_number} -v minor_idx=$i \
        '{ for (i=1; i<=NF; i++) $i=(a[i]+=$i) }END{ print cancer, major_idx, minor_idx, $0 }' \
      >> permutation_results.tmp.tsv

      # Update final output file with results from this permutation
      cat permutation_results.tmp.tsv >> ~{outfile_name}

    done
  >>>
  
  runtime {
    cpu: n_cpu
    memory: select_first([mem_gb, 1.75 * n_cpu]) + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }

  output {
    File results_tsv = outfile_name
    Array[String] cancers = read_lines("cancers.list")
  }
}


task ComparePermutedAndEmpirical {
  input {
    File observed_overlaps_tsv
    File permuted_overlaps_tsv
    Array[String] cancers

    File analysis_script
    String output_prefix

    String docker
    Int n_cpu = 4
    Float? mem_gb
    Int disk_gb = 30
  }

  command <<<
    set -eu -o pipefail

    mkdir "~{output_prefix}_results"

    # First count the number of overlaps in the empirically observed data
    while read cancer; do

      # Write one line per cancer
      for wrapper in 1; do
        echo $cancer

        # Get count of convergent gene pairs irrespective of germline|somatic context
        awk -v FS="\t" -v cancer="$cancer" \
          '{ if ($1==cancer) print }' \
          ~{observed_overlaps_tsv} | wc -l
        for criteria in same_gene ligand_receptor known_ppi protein_complex; do
          awk -v FS="\t" -v cancer="$cancer" -v criteria=$criteria \
            '{ if ($1==cancer && $NF ~ criteria) print }' \
            ~{observed_overlaps_tsv} | wc -l
        done

        # Enumerate results by germline & somatic contexts
        for germ_context in coding noncoding; do
          for som_context in coding noncoding; do
            awk -v FS="\t" -v cancer="$cancer" -v gc=$germ_context -v sc=$som_context \
              '{ if ($1==cancer && $3==gc && $5==sc) print }' \
              ~{observed_overlaps_tsv} | wc -l
            for criteria in same_gene ligand_receptor known_ppi protein_complex; do
              awk -v FS="\t" -v cancer="$cancer" -v gc=$germ_context \
                -v sc=$som_context -v criteria=$criteria \
                '{ if ($1==cancer && $3==gc && $5==sc && $NF ~ criteria) print }' \
                ~{observed_overlaps_tsv} | wc -l
            done
          done
        done

      done | paste -s

    done < ~{write_lines(cancers)} \
    > "~{output_prefix}_results/~{output_prefix}.observed_counts.tsv"

    # Add one line for all cancers
    cut -f2- "~{output_prefix}_results/~{output_prefix}.observed_counts.tsv" \
    | awk -v OFS="\t" -v cancer="all" \
      '{ for (i=1; i<=NF; i++) $i=(a[i]+=$i) }END{ print cancer, $0 }' \
    >> "~{output_prefix}_results/~{output_prefix}.observed_counts.tsv"

    # Next, compare the observed values to expected values under a permuted ull
    Rscript ~{analysis_script} \
      "~{output_prefix}_results/~{output_prefix}.observed_counts.tsv" \
      ~{permuted_overlaps_tsv} \
      "~{output_prefix}_results/~{output_prefix}"

    # Finally, create a compressed tarball of the results
    tar -czvf "~{output_prefix}_results.tar.gz" "~{output_prefix}_results"
  >>>
  
  runtime {
    cpu: n_cpu
    memory: select_first([mem_gb, 1.75 * n_cpu]) + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }

  output {
    File results_tarball = "~{output_prefix}_results.tar.gz"
  }
}


task BundleOutputs {
  input {
    Array[File] tarballs
    String output_prefix
    String docker
  }

  Int disk_gb = ceil(3 * size(tarballs, "GB")) + 20

  command <<<
    set -eu -o pipefail

    mkdir ~{output_prefix}_results

    while read tarball; do
      tar -xzvf $tarball -C ~{output_prefix}_results/
    done < ~{write_lines(tarballs)}

    tar -czvf ~{output_prefix}_results.tar.gz ~{output_prefix}_results
  >>>
  
  runtime {
    memory: "1.75 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }

  output {
    File results_tarball = "~{output_prefix}_results.tar.gz"
  }
}

