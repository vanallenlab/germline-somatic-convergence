# Germline:Somatic Convergence Project
# Copyright (c) 2023 Ryan Collins and the Van Allen Lab @ Dana-Farber Cancer Institute
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
    File gene_chrom_map_tsv
    File expression_quantiles


    # Reference files for assessing overlap
    File cellchat_csv
    File ppi_tsv
    File complexes_tsv

    # Certain scripts supplied as inputs to avoid needing to build a special
    # Docker image specifically for this analysis
    File shuffle_script
    File find_pairs_script
    File analysis_script

    # Parallelization control
    Int n_shards = 400
    Int n_perm_per_shard = 25
    Int n_cpu_per_shard = 2
    Float? mem_gb_per_shard
    Int disk_gb_per_shard = 20

    # Analysis control
    Int analysis_n_cpu = 4
    Float? analysis_mem_gb
    Int analysis_disk_gb = 30
    
    String docker = "vanallenlab/rlctools:latest"
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
        gene_chrom_map_tsv = gene_chrom_map_tsv,
        cellchat_csv = cellchat_csv,
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
        gene_chrom_map_tsv = gene_chrom_map_tsv,
        cellchat_csv = cellchat_csv,
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
        gene_chrom_map_tsv = gene_chrom_map_tsv,
        expression_quantiles = expression_quantiles,
        expression_quantile_gene_counts_tsv = expression_quantile_gene_counts_tsv,
        cellchat_csv = cellchat_csv,
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
        gene_chrom_map_tsv = gene_chrom_map_tsv,
        expression_quantiles = expression_quantiles,
        expression_quantile_gene_counts_tsv = expression_quantile_gene_counts_tsv,
        cellchat_csv = cellchat_csv,
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
    File gene_chrom_map_tsv
    File? expression_quantiles
    File? expression_quantile_gene_counts_tsv

    File cellchat_csv
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
    Int n_preemptible = 1
  }

  String outfile_name = "perm_results." + shard_number + ".tsv"

  command <<<
    set -eu -o pipefail

    # Get list of cancer types
    cut -f1 ~{strata_gene_counts_tsv} | sort -V | uniq > cancers.list

    # Run multiple permutations in each task
    for i in $( seq 1 ~{n_perm_per_shard} ); do

      date
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
      while read cancer origin context contig n_genes; do

        # Subset eligible genes to chromosome of interest
        awk -v chrom="$contig" '{ if ($2==chrom) print $1 }' ~{gene_chrom_map_tsv} \
        | fgrep -xf ~{eligible_gene_symbols} \
        > "elig.$contig.genes.list"

        # Define strata-specific weights
        case "${origin}_${context}" in
          germline_coding_gwas)
            weights=~{coding_gwas_weights}
            # For coding GWAS hits, by design we want to exclude genes that were
            # already sampled from the coding germline COSMIC set
            fgrep \
              -xvf gene_lists/germline_coding_cosmic/$cancer.germline.coding_cosmic.genes.list \
              "elig.$contig.genes.list" \
            > tmp.coding_gwas.elig_genes.list
            elig_list=tmp.coding_gwas.elig_genes.list
            ;;
          germline_noncoding)
            weights=~{noncoding_gwas_weights}
            elig_list="elig.$contig.genes.list"
            ;;
          somatic_noncoding)
            weights=~{somatic_noncoding_weights}
            elig_list="elig.$contig.genes.list"
            ;;
          *)
            weights=~{coding_nonsyn_weights}
            elig_list="elig.$contig.genes.list"
            ;;
        esac

        # Set seed string
        # Note: we take the first and last characters of each cancer, origin, 
        # and context because of total seed string length limits in our R script
        # The only exception is for negative control phenotypes, where we take 
        # the first letter of the cancer and of the non-cancer phenotype to avoid
        # a collision between mycardial_infarction and atrial_fibrilation
        sp1="~{shard_number}"
        if [ $( echo $cancer | fgrep "_" | wc -l ) -gt 0 ]; then
          sp2="$( echo $cancer | sed 's/_/\n/g' | cut -c1 | sed -n '1,2p' | paste -s -d "" )"
        else
          sp2="$( echo $cancer | cut -c1 )${cancer: -1}"
        fi
        sp3="$( echo $origin | cut -c1 )${origin: -1}"
        sp4="$( echo $context | cut -c1 )${context: -1}"
        sp5="$i"
        seed="${sp1}${sp2}${sp3}${sp4}${sp5}"

        if [ ~{defined(expression_quantiles)} == "true" ] && \
           [ ~{defined(expression_quantile_gene_counts_tsv)} == "true" ]; then
          # Shuffle & sample each expression quantile separately, if optioned
          while read q n; do
            if [ $n -eq 0 ]; then
              touch gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
              continue
            fi
            # Define list of genes in this expression quantile in this tissue
            awk -v FS="\t" -v cancer=$cancer -v q=$q \
              '{ if ($2==q && $3==cancer) print $1 }' \
              ~{expression_quantiles} \
            | fgrep -xf $elig_list \
            | sort -V > $cancer.gex_q$q.genes.list
            # Shuffle genes, restrict to eligible gene symbols, and sample desired number
            qseed="${seed}${q}"
            Rscript ~{shuffle_script} \
              --tsv-in $weights \
              --seed "$qseed" \
              --eligible-genes $cancer.gex_q$q.genes.list \
              --no-zero-weights \
              --outfile tmp.shuffled.genes.list
            # Note: need to break this up to avoid Rcript throwing SIGPIPE error
            head -n $n tmp.shuffled.genes.list \
            >> gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          done < <( awk -v FS="\t" -v OFS="\t" -v cancer="$cancer" -v origin="$origin" \
                        -v context="$context" -v contig="$contig" \
                        '{ if ($1==cancer && $3==origin && $4==context && $5==contig) print $2, $6 }' \
                        ~{expression_quantile_gene_counts_tsv} )
        else
          if [ $n_genes -gt 0 ]; then
            # Shuffle genes, restrict to eligible gene symbols, and sample desired number
            Rscript ~{shuffle_script} \
              --tsv-in $weights \
              --seed "$seed" \
              --eligible-genes $elig_list \
              --no-zero-weights \
              --outfile tmp.shuffled.genes.list
            # Note: need to break this up to avoid Rcript throwing SIGPIPE error
            head -n $n_genes tmp.shuffled.genes.list \
            >> gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          else
            touch gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
          fi
        fi

      done < ~{strata_gene_counts_tsv}

      # Consistency check: ensure number of genes sampled matches expected number
      while read cancer origin context; do
        n_expected=$( awk \
                        -v FS="\t" -v cancer="$cancer" \
                        -v origin="$origin" -v context="$context" \
                        '{ if ($1==cancer && $2==origin && $3==context) sum+=$5 }END{ print sum }' \
                        ~{strata_gene_counts_tsv} )
        n_sampled=$( cat gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list | wc -l )
        if [ $n_expected -ne $n_sampled ]; then
          echo "Error: number of sampled genes ($n_sampled) does not match number of expected genes ($n_expected)"
          echo "       for permutation $i for $cancer $origin $context. Exiting."
          exit 1
        fi
      done < <( cut -f1-3 ~{strata_gene_counts_tsv} \
                | sort -Vk1,1 -k2,2V -k3,3V | uniq )


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
              --cellchat-db ~{cellchat_csv} \
              --ppi-db ~{ppi_tsv} \
              --protein-complexes ~{complexes_tsv} \
              --report-counts
          done
        done \
        | paste -s \
        | awk -v OFS="\t" -v cancer=$cancer -v major_idx=~{shard_number} -v minor_idx=$i \
          '{ print cancer, major_idx, minor_idx, $0 }' \
        >> permutation_results.tmp.tsv

      done < cancers.list

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
    preemptible: n_preemptible
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

        # Enumerate results by germline & somatic contexts
        for germ_context in coding noncoding; do
          for som_context in coding noncoding; do
            for tier in $( seq 1 4 ); do
              awk -v FS="\t" -v cancer="$cancer" -v gc=$germ_context \
                -v sc=$som_context -v tier=$tier \
                '{ if ($1==cancer && $3==gc && $5==sc && $NF==tier) print }' \
                ~{observed_overlaps_tsv} | wc -l
            done
          done
        done

      done | paste -s

    done < ~{write_lines(cancers)} \
    > "~{output_prefix}_results/~{output_prefix}.observed_counts.tsv"

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
    preemptible: 1
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

