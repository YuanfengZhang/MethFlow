configfile: "config/runtime_config.yaml"
from textwrap import dedent
"""
! whole genome index and rrbs index cannot be generated in same dir.
BatMeth2 index -g ./BL.fa
BatMeth2 index_rrbs -g ./BL.fa
BatMeth2 align \
    -1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    -2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -g /hot_warm_data/ref/quartet/DNA/custom_genome/BL/BatMeth2/BL.fa \
    -p 16 \
    -O /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align \
    -o batmeth2
"""


def retrieve_batmeth2_ref(wildcards, config):
    match wildcards.BaseName.split('_')[0][: 2]:
        case 'RR':
            return config["ref"]["batmeth2"]["RRBS"][wildcards.BaseName.split('_')[1]]
        case _:
            return config["ref"]["batmeth2"]["default"][wildcards.BaseName.split('_')[1]]


rule batmeth2_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/batmeth2/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/batmeth2/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/batmeth2/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: retrieve_batmeth2_ref(wildcards,
                                                               config),
        extra_params = config["batmeth2"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)
        
        resources/BatMeth2/bin/BatMeth2 align \\
            -1 {input.r1} -2 {input.r2} \\
            -g {params.ref} -p {threads} \\
            -O result/{wildcards.BaseName}/{wildcards.AlignParentDir}/batmeth2 \\
            -o {wildcards.BaseName}
        
        samtools addreplacerg \\
            result/{wildcards.BaseName}/{wildcards.AlignParentDir}/batmeth2/{wildcards.BaseName}.sort.bam \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
