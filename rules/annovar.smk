### Convert vcf to avinput for annovar annotation ###

rule vcf_avinput:
    input:
        vcf_file = rules.annotate_mutect2.output['annotated_vcf'],            
    output:
        avinput = "result/030_annovar/avinput/{sample_tumor}.avinput",
    log:
        "log/030_annovar/avinput/{sample_tumor}.log"
    benchmark:
        'benchmarks/030_annovar/avinput/{sample_tumor}.benchmark'
    shell:
        """
        (convert2annovar.pl -format vcf4 -allsample -withfreq {input.vcf_file} -outfile {output.avinput}) &> {log}
        """

### Annotate with Annovar ###

rule annovar:
    input:
        avinput = rules.vcf_avinput.output['avinput']
    output:
        annovar_annotated = "result/030_annovar/{sample_tumor}.hg38_multianno.csv"
    params:
        humandb = config['resources']['humandb'],
        outfile = "result/030_annovar/{sample_tumor}"
    resources:
        num_threads = config['threads']['annovar'],
    log:
        "log/030_annovar/{sample_tumor}.log"
    benchmark:
        'benchmarks/030_annovar/{sample_tumor}.benchmark'
    shell:
        """
        (table_annovar.pl {input.avinput} {params.humandb} --buildver hg38 --outfile {params.outfile} --thread {resources.num_threads} --remove --polish --csvout --protocol refGene,ljb26_all,clinvar_20221231,revel,gnomad312_genome --operation g,f,f,f,f) &> {log}
        """

rule hgvs:
    input:
        avinput = rules.vcf_avinput.output['avinput']
    output:
        annovar_hgvs = "result/030_annovar/hgvs/{sample_tumor}.variant_function"
    params:
        humandb = config['resources']['humandb'],
        outfile = "result/030_annovar/hgvs/{sample_tumor}"
    log:
        "log/030_annovar/hgvs/{sample_tumor}.log"
    benchmark:
        'benchmarks/030_annovar/hgvs/{sample_tumor}.benchmark'
    shell:
        """
        (annotate_variation.pl -out {params.outfile} --hgvs -build hg38 {input.avinput} {params.humandb}) &> {log}
        """
