SAMPLES, = glob_wildcards("01.clean.fq/{sample}_1.fq.gz")

print("Total sample size: ",len(SAMPLES))


rule all:
    input:
        expand("01.fastp.report.summary/{sample}.csv",sample=SAMPLES)


rule fastpReportStat:
    input:
        "01.fastp.report/{sample}.json"
    output:
        "01.fastp.report.summary/{sample}.csv"
    threads:
        2
    resources:
        mem_mb = 4000
    params:
        "{sample}"
    script:
        "fastp_report_stat.R"