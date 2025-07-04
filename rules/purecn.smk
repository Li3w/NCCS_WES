### Calculate and GC-normalize coverage from a tumor BAM file ###

rule normalize_tumor:
    input:
        tumor_bam_file = "samples/tumor/{sample_tumor}.recal.bam",           
        intervals = config['resources']['intervals'],  
    output:
        tumor_GC_normalize = "result/020_normalized_purecn/tumor/{sample_tumor}.recal_coverage_loess.txt.gz",
    priority:
        1
    params:
        purecn = config['paths']['purecn'],
        out_dir_tumor_normalize = config['paths']['normalize_purecn_tumor'],
    resources:
        mem_gb = config['java_params']['custom']['mem']['mutect'] * 1024,
        cpus = config['java_params']['custom']['ncpu']['mutect']
    log:
        "log/020_normalize_purecn/{sample_tumor}.log"
    benchmark:
        'benchmarks/020_normalize_purecn/{sample_tumor}.benchmark'
    conda:
        "envs/purecn.yaml"
    shell:
        """
        (Rscript {params.purecn}/Coverage.R --out-dir {params.out_dir_tumor_normalize} --bam {input.tumor_bam_file} --intervals {input.intervals}) &> {log}
        """

rule normalize_normal:
    input:
        normal_bam_file = "samples/normal/{sample_normal}.recal.bam",
        intervals = config['resources']['intervals'],
    output:
        normal_GC_normalize = "result/020_normalized_purecn/normal/{sample_normal}.recal_coverage_loess.txt.gz"
    priority:
        1
    params:
        purecn = config['paths']['purecn'],
        out_dir_normal_normalize = config['paths']['normalize_purecn_normal'],
    resources:
        mem_gb = config['java_params']['custom']['mem']['mutect'] * 1024,
        cpus = config['java_params']['custom']['ncpu']['mutect']
    log:
        "log/020_normalize_purecn/{sample_normal}.log"
    benchmark:
        'benchmarks/020_normalize_purecn/{sample_normal}.benchmark'
    conda:
        "envs/purecn.yaml"
    shell:
        """
        (Rscript {params.purecn}/Coverage.R --out-dir {params.out_dir_normal_normalize} --bam {input.normal_bam_file} --intervals {input.intervals}) &> {log}
        """

### With a matched normal tumor sample ###
rule purecn:
    input:
        normal_interval = lambda wildcards: f"result/020_normalized_purecn/normal/{lookup_matched_normal(wildcards.sample_tumor, tumor_normal_map)}.recal_coverage_loess.txt.gz",
        tumor_interval = rules.normalize_tumor.output['tumor_GC_normalize'],
        vcf = rules.annotate_mutect2.output['annotated_vcf'],
        intervals = config['resources']['intervals'],
        stats = "result/011_mutect2_variant_filtered/stats/{sample_tumor}_mutect2_filtering.stats"
    output:
        purecn_result = "result/021_purecn/{sample_tumor}.csv"
    params:
        purecn = config['paths']['purecn'],
        tumor_id = lambda wildcards: wildcards.sample_tumor,
        out_dir_purecn = config['paths']['purecn_result']
    resources:
        mem_gb = config['java_params']['custom']['mem']['purecn'] * 1024,
        cpus = config['java_params']['custom']['ncpu']['purecn']
    log:
        "log/021_purecn/{sample_tumor}.log"
    benchmark:
        'benchmarks/021_purecn/{sample_tumor}.benchmark'
    conda:
        "envs/purecn.yaml"
    shell:
        """
        (Rscript {params.purecn}/PureCN.R --out {params.out_dir_purecn} --tumor {input.tumor_interval} --normal {input.normal_interval} --sampleid {params.tumor_id} --vcf {input.vcf} --intervals {input.intervals} --fun-segmentation PSCBS --stats-file {input.stats} --min-base-quality 20 --bootstrapn 500 --genome hg38 --force --post-optimize --seed 42) &> {log}
        """
