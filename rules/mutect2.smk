def get_tumor_normal_pairs(tsv_file):
  """
  This function reads a TSV file containing sample pairings (tumor and normal)
  and returns a dictionary mapping tumor names to their corresponding normal names.

  Args:
      tsv_file: The path to the TSV file containing sample pairings (string).

  Returns:
      A dictionary mapping tumor names to their corresponding normal names (dict).
  """

  tumor_normal_pairs = {}
  with open(tsv_file, 'r') as f:
    # Read the header line (optional, store column names if needed)
    header = next(f).strip().split('\t')
    tumor_col_index = None
    normal_col_index = None

    # Find column indexes for tumor and normal (adjust column names if needed)
    for i, col_name in enumerate(header):
      if col_name.lower() == "tumor":
        tumor_col_index = i
      elif col_name.lower() == "normal":
        normal_col_index = i

    if tumor_col_index is None or normal_col_index is None:
      raise ValueError("TSV file must contain 'tumor' and 'normal' columns")

    # Read remaining lines and build the dictionary
    for line in f:
      tumor, normal = line.strip().split('\t')
      tumor_normal_pairs[tumor] = normal
  return tumor_normal_pairs

# Call this function once (outside the rule) and store the result
tumor_normal_map = get_tumor_normal_pairs("sample_list.tsv")

def lookup_matched_normal(tumor_name, tumor_normal_map):
  """
  This function retrieves the matched normal name for a tumor sample
  from a pre-built dictionary.

  Args:
      tumor_name: The name of the tumor sample (string).
      tumor_normal_map: A dictionary mapping tumor names to normal names (dict).

  Returns:
      The matched normal name (string) or None if not found.
  """
  return tumor_normal_map.get(tumor_name, None)



### VARIANT CALLING ###

rule mutect2:
    input:
        normal_bam = lambda wildcards: f"samples/normal/{lookup_matched_normal(wildcards.sample_tumor, tumor_normal_map)}.recal.bam",
        tumor_bam = "samples/tumor/{sample_tumor}.recal.bam",
        ref = config['resources']['reference'], # fasta file
        intervals = config['resources']['bed'],  # if applicable
        gnomad = config['resources']['known_variants']['gnomad'],
    output:
        vcf = "result/010_mutect2_variant_raw/{sample_tumor}_mutect2.vcf",
        f1r2 = "result/012_f1r2/{sample_tumor}_f1r2.tar.gz",
        bam_output = "result/013_mutect2_bam/{sample_tumor}.bam",
    priority:
        7
    params:
        tumor_lod = config['params']['mutect2']['lod'],
        genotype_germline_sites = config['params']['mutect2']['genotype_germline_sites'],
        genotype_pon_sites = config['params']['mutect2']['genotype_pon_sites'],
        tumor_id = lambda wildcards: wildcards.sample_tumor,  
        normal_id = lambda wildcards: {lookup_matched_normal(wildcards.sample_tumor, tumor_normal_map)},
    resources:
        mem_mb = config['java_params']['custom']['mem']['mutect'] * 1024,
        cpus = config['java_params']['custom']['ncpu']['mutect'],
        parallel_mutect2 = 1
    log:
        "log/mutect2/{sample_tumor}.log"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        (gatk --java-options "-Xmx{resources.mem_mb}M" Mutect2 -R {input.ref} -I {input.tumor_bam} -tumor {params.tumor_id} -I {input.normal_bam} -normal {params.normal_id} --germline-resource {input.gnomad} --intervals {input.intervals} --f1r2-tar-gz {output.f1r2} --bam-output {output.bam_output} --tumor-lod-to-emit {params.tumor_lod} --genotype-germline-sites {params.genotype_germline_sites} --genotype-pon-sites {params.genotype_pon_sites} -O {output.vcf}) &> {log}
        """

### Make statistics and model about fw/rv reads ###

rule learn_read_orientation_model:
    input:
        f1r2 = "result/012_f1r2/{sample_tumor}_f1r2.tar.gz"
    output:
        artifact_prior = "result/014_artifact_prior/{sample_tumor}_artifact_prior.tar.gz"
    priority:
        6
    params:
        java_opts = config['java_params']['custom']['ncpu']['learn_read_orientation'],
        custom_mem = config['java_params']['custom']['mem']['learn_read_orientation'],
    resources:
        cpus = config['java_params']['custom']['ncpu']['learn_read_orientation'],
        mem_mb = config['java_params']['custom']['mem']['learn_read_orientation'] * 1024
    log:
        'log/read_orientation_model/{sample_tumor}.log'
    benchmark:
        'benchmarks/learn_read_orientation_model/{sample_tumor}.benchmark'
    conda:
        "envs/gatk.yaml"
    shell:
        """
        (gatk --java-options "-Xmx{resources.mem_mb}M" LearnReadOrientationModel -I {input.f1r2} -O {output.artifact_prior}) &> {log}
        """

rule filter_mutect_calls:
    input:
        unfiltered_vcf = rules.mutect2.output['vcf'],
        ref = config['resources']['reference'],
        ob_priors = rules.learn_read_orientation_model.output['artifact_prior'],
    output:
        filtered_vcf = "result/011_mutect2_variant_filtered/{sample_tumor}_mutect2.filtered.vcf",
        filtering_stats = "result/011_mutect2_variant_filtered/stats/{sample_tumor}_mutect2_filtering.stats",
    priority:
        5
    params:
        java_opts = config['java_params']['custom']['ncpu']['filter_mutect_calls'],
        custom_mem = config['java_params']['custom']['mem']['filter_mutect_calls'],
    resources:
        cpus = config['java_params']['custom']['ncpu']['filter_mutect_calls'],
        mem_mb = config['java_params']['custom']['mem']['filter_mutect_calls'] * 1024
    log:
        'log/011_filter_mutect_calls/{sample_tumor}.log'
    benchmark:
        'benchmarks/011_filter_mutect_calls/{sample_tumor}.benchmark'
    conda:
        "envs/gatk.yaml"
    shell:
        """
        (gatk FilterMutectCalls --java-options "-Xmx{resources.mem_mb}M" -V {input.unfiltered_vcf} -R {input.ref} --ob-priors {input.ob_priors} -O {output.filtered_vcf}  --filtering-stats {output.filtering_stats}) &> {log}
        """

### Annotate with DB flag ###

rule zip:
    input:
        vcf_filtered = rules.filter_mutect_calls.output['filtered_vcf']
    output:
        zip_vcf = "result/011_mutect2_variant_filtered/zip/{sample_tumor}_mutect2_filtered.vcf.gz"
    priority:
        4
    shell:
        """
        bgzip -c {input.vcf_filtered} > {output.zip_vcf} 
        """

rule index:
    input:
        zip_vcf = rules.zip.output['zip_vcf']
    output:
        index_vcf = "result/011_mutect2_variant_filtered/zip/{sample_tumor}_mutect2_filtered.vcf.gz.tbi"
    priority:
        3
    shell:
        """
        tabix -p vcf {input.zip_vcf}  
        """

rule annotate_mutect2:
    input:
        vcf_zip_filtered = rules.zip.output['zip_vcf'],
        db_snps = config['resources']['known_variants']['dbsnp'],
        index_vcf = "result/011_mutect2_variant_filtered/zip/{sample_tumor}_mutect2_filtered.vcf.gz.tbi",
    output:
        annotated_vcf = "result/015_mutect2_annotate/{sample_tumor}_mutect2_filtered_annotated.vcf.gz",
    priority:
        1
    params:
        java_opts = config['java_params']['custom']['ncpu']['annotate_mutect'],
        custom_mem = config['java_params']['custom']['mem']['annotate_mutect'],
    resources:
        cpus = config['java_params']['custom']['ncpu']['annotate_mutect'],
        mem_mb = config['java_params']['custom']['mem']['annotate_mutect'] * 1024
    log:
        'log/015_mutect2_annotate/{sample_tumor}.log'
    benchmark:
        'benchmarks/015_mutect2_annotate/{sample_tumor}.benchmark'
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        (bcftools annotate --annotation {input.db_snps} --columns ID,INFO --output {output.annotated_vcf} {input.vcf_zip_filtered}) &> {log}
        """


