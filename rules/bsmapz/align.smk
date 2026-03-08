configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
! The length of lines in the .fa is limited to 70 else malloc error pops.
! seqtk seq -l 70 BL.fa > BL/bsmapz/BL.fa
bsmapz \\
    -a /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    -b /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -d /hot_warm_data/ref/quartet/DNA/custom_genome/BL/bsmapz/BL.fa \
    -o /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/bsmapz.bam \
    -p 16

! It's has to be done in the input dir since long file paths cause the following error:
!    **.bam.tmp.bam*** buffer overflow detected ***: terminated
! see https://github.com/zyndagj/BSMAPz/issues/19
"""


def retrieve_method_params(wildcards, three: bool = False):
    to_return = ""
    if wildcards.BaseName.split('_')[0][: 2] == "RR":
        to_return += "-s 12"

    if three:
        to_return += " -3"

    return to_return


rule bsmapz_align:
    input:
        r1            = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam           = "result/{BaseName}/{AlignParentDir}/bsmapz/{BaseName}.bam",
        bai           = "result/{BaseName}/{AlignParentDir}/bsmapz/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bsmapz/{BaseName}.align.benchmark"
    params:
        ref           = lambda wildcards: config["ref"]["bsmapz"][wildcards.BaseName.split('_')[1]],
        method_params = lambda wildcards: retrieve_method_params(wildcards),
        extra_params  = config["bsmapz"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        cd result/{wildcards.BaseName}/{wildcards.AlignParentDir}
        
        ln -sf {params.ref} $(pwd)/ref.fa

        bsmapz \\
            -a {wildcards.BaseName}.R1.fq.gz \\
            -b {wildcards.BaseName}.R2.fq.gz -d ref.fa \\
            -p {threads} {params.extra_params} |\\
        mbuffer -m 4G -q |\\
        samtools sort -@ {threads} - |\\
        mbuffer -m 4G -q |\\
        samtools addreplacerg \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o bsmapz/{wildcards.BaseName}.bam -

        samtools index \\
            -@ {threads} \\
            bsmapz/{wildcards.BaseName}.bam || echo "suppress non-zero exit"

        rm ref.fa
        """)
