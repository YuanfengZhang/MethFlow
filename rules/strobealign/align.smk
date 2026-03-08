configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
strobealign --create-index BL.fa -t 16 -r 150
strobealign \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/strobealign/BL.fa \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    --use-index --mcs -t 24 \
    --rg-id ID:1 --rg SM:sample1 --rg PL:ILLUMINA --rg PU:unit1 --rg LB:lib1 |\
samtools sort \
    --output-fmt bam,level=9 \
    -@ 24 >/mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/strobealign.bam
"""


rule strobealign:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/strobealign/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/strobealign/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/strobealign/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["strobealign"][wildcards.BaseName.split('_')[1]],
        extra_params = config["strobealign"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        strobealign \\
            {params.ref} {input.r1} {input.r2} \\
            --use-index -t {threads} {params.extra_params} \\
            --rg-id ID:{wildcards.BaseName} \\
            --rg SM:${{SAMPLE}} \\
            --rg PL:${{PLATFORM}} \\
            --rg LB:${{LIB}} |\\
        mbuffer -m 4G -q |\\
        samtools sort - \\
            --output-fmt bam,level=9 \\
            -@ {threads} > {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
