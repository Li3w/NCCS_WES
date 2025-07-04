### Calculate Microsatalite Instability with MSISensor-pro ###

rule msi:
    input:
        normal_bam = lambda wildcards: f"samples/normal/{lookup_matched_normal(wildcards.sample_tumor, tumor_normal_map)}.recal.bam",
        tumor_bam = "samples/tumor/{sample_tumor}.recal.bam",
        msi_ref = config['resources']['MSI']
    output:
        msi = "result/060_msi/{sample_tumor}"
    log:
        "log/060_msi/{sample_tumor}.log"
    benchmark:
        'benchmarks/060_msi/{sample_tumor}.benchmark'
    conda:
        "envs/msi.yaml"
    shell:
        """
        (msisensor-pro msi -d {input.msi_ref} -n {input.normal_bam} -t {input.tumor_bam} -o {output}) &> {log}
        """
