### Convert vcf to maf for oncokb annotation ###

rule gunzip:
    input:
        vcf_file = rules.annotate_mutect2.output['annotated_vcf']
    output:
        vcf_file_unzip = "result/015_mutect2_annotate/{sample_tumor}_mutect2_filtered_annotated.vcf"
    shell:
        """
        gunzip -c {input.vcf_file} > {output.vcf_file_unzip}
        """

rule vcf_maf:
    input:
        vcf_file_unzip = rules.gunzip.output['vcf_file_unzip'],
        ref = config['resources']['reference'], #fasta file           
    output:
        maf = "result/040_oncokb/maf/{sample_tumor}.maf",
    params:
        vcf2maf = config['paths']['vcf2maf'],
        tumor_id = lambda wildcards: wildcards.sample_tumor,  
        normal_id = lambda wildcards: {lookup_matched_normal(wildcards.sample_tumor, tumor_normal_map)},
        vep_path = config['paths']['vep'],
        ncbi_build = "GRCh38"
    resources:
        cpus = config['threads']['vcf2maf'],
        mem_mb = config['resources']['mem']['vcf2maf'] * 1024,
    log:
        "log/040_oncokb/maf/{sample_tumor}.log"
    benchmark:
        'benchmarks/040_oncokb/maf/{sample_tumor}.benchmark'
    conda:
        "envs/oncokb.yaml"
    shell:
        """
        (perl {params.vcf2maf}/vcf2maf.pl --input-vcf {input.vcf_file_unzip} --output-maf {output.maf} --ref-fasta {input.ref} --vep-forks {resources.cpus} --inhibit-vep --tumor-id {params.tumor_id} --normal-id {params.normal_id} --ncbi-build {params.ncbi_build}) &> {log}
        """

### Annotate with oncokb ###

rule oncokb:
    input:
        maf = rules.vcf_maf.output['maf']
    output:
        oncokb_annotated = "result/040_oncokb/{sample_tumor}_oncokb.maf"
    params:
        oncokb_path = config['paths']['oncokb'],
        token = config['params']['oncokb']['token'],
    #resources:
        #num_threads = config['threads']['annovar'],
    log:
        "log/040_oncokb/{sample_tumor}.log"
    benchmark:
        'benchmarks/040_oncokb/{sample_tumor}.benchmark'
    conda:
        "envs/oncokb.yaml"
    shell:
        """
        (python {params.oncokb_path}/MafAnnotator.py -i {input.maf} -o {output.oncokb_annotated} -q genomic_change -b {params.token}) &> {log}
        """

