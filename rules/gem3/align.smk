configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
./configure && make HAVE_CUDA=0 all
gem-indexer -i BL.fa -o BL_bisulfite -b -t 64 -v
gem-indexer -i BL.fa -o BL -t 64 -v
gem-mapper \
    -I /hot_warm_data/ref/quartet/DNA/custom_genome/BL/gem3/BL_bisulfite.gem \
    -1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.R1.fq.gz \
    -2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.R2.fq.gz \
    -r "@RG\tID:1\tSM:sample1\tPL:ILLUMINA\tPU:unit1\tLB:lib1" \
    -p -z -t 64 |\
samtools sort \
    --output-fmt bam,level=9 \
    -@ 64 \
    -o /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/gem3_bisulfite.bam
"""


def retrieve_gem3_ref(wildcards, config):
    if wildcards.BaseName.split('_')[0][: 2] == 'PS':
        return config["ref"]["gem3"]["original"][wildcards.BaseName.split('_')[1]]
    else:
        return config["ref"]["gem3"]["bisulfite"][wildcards.BaseName.split('_')[1]]


rule gem3_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/gem3/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/gem3/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/gem3/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: retrieve_gem3_ref(wildcards,
                                                           config),
        extra_params = config["gem3"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        resources/gem3-mapper/bin/gem-mapper \\
            -I {params.ref} -1 {input.r1} -2 {input.r2} \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            -p -z {params.extra_params} \\
            -t {threads} |\\
        samtools sort \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
