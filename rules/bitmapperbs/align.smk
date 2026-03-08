configfile: "config/runtime_config.yaml"
from textwrap import dedent
"""
bitmapperBS --index BL.fa
bitmapperBS search \\
    -t 8 --sensitive \\
    --search /hot_warm_data/ref/quartet/DNA/custom_genome/BL/BitMapperBS/BL.fa \\
    --seq1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \\
    --seq2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz |\\
    samtools sort -O bam \\
    -o /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/biitmapperbs.bam
"""


rule bitmapperbs_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bitmapperbs/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bitmapperbs/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bitmapperbs/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bitmapperbs"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bitmapperbs"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        bitmapperBS search \\
            -t {threads} {params.extra_params} \\
            --search {params.ref} \\
            --seq1 {input.r1} --seq2 {input.r2} |\\
        mbuffer -m 4G -q |\\
        samtools sort \\
            -@ {threads} |\\
        samtools addreplacerg - \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
