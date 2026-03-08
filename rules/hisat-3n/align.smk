configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout -b hisat-3n origin/hisat-3n
make
time hisat-3n-build \
    ./BL.fa ./BL_c2t \
    --base-change C,T \
    -p 24
time hisat-3n-build \
    ./BL.fa ./BL_t2c \
    --base-change T,C \
    -p 24 \
    --repeat-index
hisat-3n \
    --index /hot_warm_data/ref/quartet/DNA/custom_genome/BL/hisat-3n/BL_c2t \
    -q \
    -1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    -2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -S /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/hisat3n.sam \
    --base-change C,T --sensitive \
    --rg-id ID:1 \
    --rg SM:sample1 \
    --rg PL:ILLUMINA \
    --rg PU:unit1 \
    --rg LB:lib1 \
    -t -p 16
"""

def retrieve_hisat3n_ref(wildcards, config):
    if wildcards.BaseName.split('_')[0][: 2] == 'PS':
        return config["ref"]["hisat-3n"][wildcards.BaseName.split('_')[1]]["t2c"]
    else:
        return config["ref"]["hisat-3n"][wildcards.BaseName.split('_')[1]]["c2t"]


def retrieve_base_change(wildcards):
    if wildcards.BaseName.split('_')[0][: 2] == 'PS':
        return "T,C"
    else:
        return "C,T"


rule hisat3n_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/hisat-3n/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/hisat-3n/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/hisat-3n/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: retrieve_hisat3n_ref(wildcards,
                                                               config),
        base_change  = lambda wildcards: retrieve_base_change(wildcards),
        extra_params = config["hisat-3n"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        resources/hisat-3n/hisat-3n \\
            --index {params.ref} \\
            -q -1 {input.r1} -2 {input.r2} \\
            --base-change {params.base_change} \\
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
