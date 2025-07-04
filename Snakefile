import pandas as pd

sample_list = pd.read_csv("sample_list.tsv", sep = "\t")

#### config 
configfile: "config.yaml"

#### load rule
include: "rules/mutect2.smk"
include: "rules/purecn.smk"
include: "rules/annovar.smk"
include: "rules/oncokb.smk"
include: "rules/pytmb.smk"
include: "rules/msi.smk"

### SOMATIC TUMOR-NORMAL VARIANT CALLING ###
rule all:
    input:
        expand("result/010_mutect2_variant_raw/{sample_tumor}_mutect2.vcf", sample_tumor = sample_list['tumor']),
        expand("result/014_artifact_prior/{sample_tumor}_artifact_prior.tar.gz", sample_tumor = sample_list['tumor']),
        expand("result/011_mutect2_variant_filtered/{sample_tumor}_mutect2.filtered.vcf", sample_tumor = sample_list['tumor']),
        expand("result/011_mutect2_variant_filtered/zip/{sample_tumor}_mutect2_filtered.vcf.gz.tbi", sample_tumor = sample_list['tumor']),
        expand("result/015_mutect2_annotate/{sample_tumor}_mutect2_filtered_annotated.vcf.gz", sample_tumor = sample_list['tumor']),
        expand("result/020_normalized_purecn/tumor/{sample_tumor}.recal_coverage_loess.txt.gz", sample_tumor = sample_list['tumor']),
        expand("result/020_normalized_purecn/normal/{sample_normal}.recal_coverage_loess.txt.gz", sample_normal = sample_list['normal']),
        expand("result/021_purecn/{sample_tumor}.csv", sample_tumor = sample_list['tumor']),
        expand("result/030_annovar/{sample_tumor}.hg38_multianno.csv", sample_tumor = sample_list['tumor']),
        expand("result/030_annovar/hgvs/{sample_tumor}.variant_function", sample_tumor = sample_list['tumor']),
        expand("result/040_oncokb/maf/{sample_tumor}.maf", sample_tumor = sample_list['tumor']),
        expand("result/040_oncokb/{sample_tumor}_oncokb.maf", sample_tumor = sample_list['tumor']),
	    expand("result/050_tmb/{sample_tumor}.txt", sample_tumor = sample_list['tumor']),
	    expand("result/060_msi/{sample_tumor}", sample_tumor = sample_list['tumor'])
