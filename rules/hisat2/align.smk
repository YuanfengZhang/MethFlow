configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
hisat2-build BL.fa BL
hisat2 -x /hot_warm_data/ref/quartet/DNA/custom_genome/BL/hisat2/BL \
    -1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    -2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -q --sensitive \
    --rg-id ID:1 \
    --rg SM:sample1 \
    --rg PL:ILLUMINA \
    --rg PU:unit1 \
    --rg LB:lib1 \
    -t -p 16
"""


rule hisat2_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/hisat2/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/hisat2/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/hisat2/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["hisat2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["hisat2"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        resources/hisat2/hisat2 \\
            -x {params.ref} -1 {input.r1} -2 {input.r2} \\
            --rg-id ID:{wildcards.BaseName} \\
            --rg SM:${{SAMPLE}} \\
            --rg PL:${{PLATFORM}}\\
            --rg LB:${{LIB}} \\
            -p {threads} {params.extra_params} |\\
        mbuffer -m 4G -q |\\
        samtools sort \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
