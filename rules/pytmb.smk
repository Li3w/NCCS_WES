### Calculate Tumor mutational burden (TMB) ###

# rule genome_size:
#     input:
#         bam = rules.mutect2.output['bam_output']
#     output:
#         genome_size = "result/050_tmb/genome_size/{sample_tumor}.thresholds.bed"
#     params:
#         bed = config['resources']['bed_sorted'],
#         gencode = config['resources']['gencode']
#     log:
#         "log/050_tmb/genome_size/{sample_tumor}.log"
#     benchmark: 
#         'benchmarks/050_tmb/genome_size/{sample_tumor}.benchmark'
#     conda:
#         "envs/tmb.yaml"  
#     shell:
#         """
#         (python TMB-master/bin/pyEffGenomeSize.py --bed {params.bed} --gtf {params.gencode} --bam {input.bam} --oprefix {wildcards.sample_tumor} --mosdepth --filterNonCoding --minCoverage 10) &> {log}
#         """

rule tmb:
    input:
        vcf = rules.annotate_mutect2.output['annotated_vcf']
    output:
        tmb = "result/050_tmb/{sample_tumor}.txt"
    params:
        vaf = config['params']['tmb']['vaf'],
        maf = config['params']['tmb']['maf'],
        mindepth = config['params']['tmb']['min_depth'],
    resources:
        dbconfig = config['resources']['TMB']['dbconfig'],
        varconfig = config['resources']['TMB']['varconfig']
    log:
        "log/050_tmb/{sample_tumor}.log"
    benchmark:
        'benchmarks/050_tmb/{sample_tumor}.benchmark'
    conda:
        "envs/tmb.yaml"
    shell:
        """
        (python TMB-master/bin/pyTMB.py -i {input.vcf} --effGenomeSize 33280000 --dbConfig {resources.dbconfig} --varConfig {resources.varconfig} --vaf {params.vaf} --maf {params.maf} --minDepth {params.mindepth} --sample {wildcards.sample_tumor} --minAltDepth 2 --filterLowQual --filterNonCoding --filterRecurrence --filterSyn --filterCancerHotspot --filterPolym --polymDb 1k,gnomad > {output.tmb}) &> {log}
        """
