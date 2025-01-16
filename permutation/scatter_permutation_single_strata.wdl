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
    File postprocess_script
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
    call PostprocessOverlaps as PostprocessUniform {
      input:
        tsv_in = PermuteUniform.results_tsv,
        postprocess_script = postprocess_script,
        docker = docker
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
    call PostprocessOverlaps as PostprocessBayesian {
      input:
        tsv_in = PermuteBayesian.results_tsv,
        postprocess_script = postprocess_script,
        docker = docker
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
    call PostprocessOverlaps as PostprocessExpression {
      input:
        tsv_in = PermuteExpression.results_tsv,
        postprocess_script = postprocess_script,
        docker = docker
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
    call PostprocessOverlaps as PostprocessComposite {
      input:
        tsv_in = PermuteComposite.results_tsv,
        postprocess_script = postprocess_script,
        docker = docker
    }
  }

  # Concatenate results from uniformly weighted permutations and compare to empirically observed results
  call Utils.ConcatTextFiles as ConcatUniform {
    input:
      shards = PostprocessUniform.postprocessed_results,
      compression_command = "gzip -c",
      output_filename = output_prefix + ".permutation_results.uniform_sampling.n" + total_n_perms + ".txt.gz",
      docker = docker
  }
  call ComparePermutedAndEmpirical as AnalyzeUniform {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatUniform.merged_file,
      cancers = PermuteUniform.cancers[0],
      analysis_script = analysis_script,
      output_prefix = output_prefix + ".uniform_sampling",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from Bayesian weighted permutations and compare to empirically observed results
  call Utils.ConcatTextFiles as ConcatBayesian {
    input:
      shards = PostprocessBayesian.postprocessed_results,
      compression_command = "gzip -c",
      output_filename = output_prefix + ".permutation_results.bayesian_sampling.n" + total_n_perms + ".txt.gz",
      docker = docker
  }
  call ComparePermutedAndEmpirical as AnalyzeBayesian {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatBayesian.merged_file,
      cancers = PermuteBayesian.cancers[0],
      analysis_script = analysis_script,
      output_prefix = output_prefix + ".bayesian_sampling",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from expression quantile-matched permutations and compare to empirically observed results
  call Utils.ConcatTextFiles as ConcatExpression {
    input:
      shards = PostprocessExpression.postprocessed_results,
      compression_command = "gzip -c",
      output_filename = output_prefix + ".permutation_results.expression_sampling.n" + total_n_perms + ".txt.gz",
      docker = docker
  }
  call ComparePermutedAndEmpirical as AnalyzeExpression {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatExpression.merged_file,
      cancers = PermuteExpression.cancers[0],
      analysis_script = analysis_script,
      output_prefix = output_prefix + ".expression_sampling",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Concatenate results from composite weighted permutations and compare to empirically observed results
  call Utils.ConcatTextFiles as ConcatComposite {
    input:
      shards = PostprocessComposite.postprocessed_results,
      compression_command = "gzip -c",
      output_filename = output_prefix + ".permutation_results.composite_sampling.n" + total_n_perms + ".txt.gz",
      docker = docker
  }
  call ComparePermutedAndEmpirical as AnalyzeComposite {
    input:
      observed_overlaps_tsv = observed_overlaps_tsv,
      permuted_overlaps_tsv = ConcatComposite.merged_file,
      cancers = PermuteComposite.cancers[0],
      analysis_script = analysis_script,
      output_prefix = output_prefix + ".composite_sampling",
      docker = docker,
      n_cpu = analysis_n_cpu,
      mem_gb = analysis_mem_gb,
      disk_gb = analysis_disk_gb
  }

  # Bundle all results as a single tarball for returning to parent workflow
  call BundleOutputs {
    input:
      tarballs = [AnalyzeUniform.results_tarball, AnalyzeBayesian.results_tarball, 
                  AnalyzeExpression.results_tarball, AnalyzeComposite.results_tarball],
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
      if [ -e perm_univ ]; then
        rm -rf perm_univ
      fi
      mkdir perm_univ
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
      # Do this up front for all unique strata before counting
      echo -e "\nNow permuting all eligible genes for all strata. Seeds:"
      time while read cancer origin context; do
    
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
        echo "$seed"

        # Define strata-specific weights
        case "${origin}_${context}" in
          germline_coding_gwas)
            weights=~{coding_gwas_weights}
            ;;
          germline_noncoding)
            weights=~{noncoding_gwas_weights}
            ;;
          somatic_noncoding)
            weights=~{somatic_noncoding_weights}
            ;;
          *)
            weights=~{coding_nonsyn_weights}
            ;;
        esac

        # Shuffle the universe of all genes, which can be subset and sampled later
        Rscript ~{shuffle_script} \
          --tsv-in "$weights" \
          --seed "$seed" \
          --eligible-genes ~{eligible_gene_symbols} \
          --no-zero-weights \
          --outfile "perm_univ/$cancer.$origin.$context.shuffled.genes.list"

        # Split shuffled universe by contig for convenience
        while read contig; do
            awk -v chrom="$contig" '{ if ($2==chrom) print $1 }' ~{gene_chrom_map_tsv} \
            | fgrep -xf - "perm_univ/$cancer.$origin.$context.shuffled.genes.list" \
            > "perm_univ/$cancer.$origin.$context.$contig.shuffled.genes.list"
        done < <( cut -f4 ~{strata_gene_counts_tsv} | sort -V | uniq )

      done < <( cut -f1-3 ~{strata_gene_counts_tsv} \
                | sort -Vk1,1 -k2,2V -k3,3V | uniq )

      # Sample gene lists from permuted universe for each strata
      echo -e "\nNow sampling from permuted gene lists"
      time while read cancer origin context contig n_genes; do

        # For coding GWAS hits, by design we want to exclude genes that were
        # already sampled from the coding germline COSMIC set
        contig_univ="perm_univ/$cancer.$origin.$context.$contig.shuffled.genes.list"
        case "${origin}_${context}" in
          germline_coding_gwas)
            fgrep \
              -xvf gene_lists/germline_coding_cosmic/$cancer.germline.coding_cosmic.genes.list \
              "$contig_univ" \
            > tmp.coding_gwas.elig_genes.list || true
            elig_list=tmp.coding_gwas.elig_genes.list
            ;;
          *)
            elig_list="$contig_univ"
            ;;
        esac

        # Sample genes from pre-shuffled universe
        target=gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list
        touch $target
        if [ $n_genes -gt 0 ]; then
          if [ ~{defined(expression_quantiles)} == "true" ] && \
             [ ~{defined(expression_quantile_gene_counts_tsv)} == "true" ]; then
            
            # Sample each expression quantile separately, if optioned
            while read q n; do
              if [ $n -gt 0 ]; then
                # Further restrict eligible genes list to this expression quantile in this tissue
                awk -v FS="\t" -v cancer=$cancer -v q=$q \
                  '{ if ($2==q && $3==cancer) print $1 }' \
                  ~{expression_quantiles} \
                | fgrep -xf - $elig_list \
                > $cancer.gex_q$q.genes.list || true
                n_elig=$( cat $cancer.gex_q$q.genes.list | wc -l )
                if [ $n_elig -gt 0 ]; then
                  head -n $n $cancer.gex_q$q.genes.list >> $target || true
                fi
              fi
            done < <( awk -v FS="\t" -v OFS="\t" -v cancer="$cancer" -v origin="$origin" \
                          -v context="$context" -v contig="$contig" \
                          '{ if ($1==cancer && $3==origin && $4==context && $5==contig) print $2, $6 }' \
                          ~{expression_quantile_gene_counts_tsv} )
          else
            head -n $n_genes $elig_list \
            >> gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list || true
          fi
        fi

      done < ~{strata_gene_counts_tsv}

      # Consistency check: ensure number of genes sampled matches expected number
      echo -e "\nChecking consistency of sampled and expected gene lists"
      time while read cancer origin context; do
        # Ensure deduplication of sampled list
        sort -V gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list \
        | uniq > gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list2
        mv \
          gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list2 \
          gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list

        # Compare sizes of sampled and expected lists
        n_expected=$( awk \
                        -v FS="\t" -v cancer="$cancer" \
                        -v origin="$origin" -v context="$context" \
                        '{ if ($1==cancer && $2==origin && $3==context) sum+=$5 }END{ print sum }' \
                        ~{strata_gene_counts_tsv} )
        n_sampled=$( cat gene_lists/${origin}_${context}/$cancer.$origin.$context.genes.list | wc -l )
        if [ $n_expected -ne $n_sampled ]; then
          echo "Error: number of sampled genes ($n_sampled) does not match number of expected genes ($n_expected) for permutation $i for $cancer $origin $context. Exiting."
          exit 1
        fi
      done < <( cut -f1-3 ~{strata_gene_counts_tsv} \
                | sort -Vk1,1 -k2,2V -k3,3V | uniq )


      # Find pairs per cancer type
      echo -e "\nFinding convergent pairs in permuted gene lists"
      time while read cancer; do

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


task PostprocessOverlaps {
  input {
    File tsv_in
    File postprocess_script

    String docker

    Int n_cpu = 2
    Float? mem_gb
    Int n_preemptible = 1
  }

  String outfile_name = basename(tsv_in, ".tsv") + ".postprocessed.tsv"
  Int disk_gb = ceil(3 * size(tsv_in, "GiB")) + 10

  command <<<
    set -eu -o pipefail

    Rscript ~{postprocess_script} ~{tsv_in} ~{outfile_name}
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
    File postprocessed_results = outfile_name
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

