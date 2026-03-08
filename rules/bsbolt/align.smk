configfile: "config/runtime_config.yaml"
from textwrap import dedent
"""
bsbolt Index -G ./BL.fa -DB .
bsbolt Align \\
    -F1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    -F2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -O /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/BS_HF_1_MX \
    -DB /hot_warm_data/ref/quartet/DNA/custom_genome/BL/BSBolt2 \
    -R "@RG ID:1\tSM:sample1\tPL:ILLUMINA\tPU:unit1\tLB:lib1" \
    -UN -OT 8 -t 8

# ! the bsbolt requires outdated libcrypto.so.1.1 with is incompatible with samtools.
# ! That's why the bsbolt align and post-align are separated.

# ! bsbolt won't work within docker image and job system, but only with host machine.
# ! Somethings is written wrong in the Externel/BWA folder.
# ! See also: https://github.com/NuttyLogic/BSBolt/issues/158
"""

rule bsbolt_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bsbolt/tmp_{BaseName}.bam"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bsbolt/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bsbolt"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bsbolt"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        target_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bsbolt"
    
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        mkdir -p ${{target_dir}}/tmp_{wildcards.BaseName}

        bsbolt Align \\
            -F1 {input.r1} -F2 {input.r2} \\
            -O ${{target_dir}}/tmp_{wildcards.BaseName} \\
            -DB {params.ref} \\
            -R "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            -UN -OT {threads} -t {threads}
        """)

rule bsbolt_post_align:
    input:
        "result/{BaseName}/{AlignParentDir}/bsbolt/tmp_{BaseName}.bam"
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bsbolt/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bsbolt/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bsbolt/{BaseName}.post_align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bsbolt"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bsbolt"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        target_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bsbolt"

        samtools sort \\
            --output-fmt bam,level=9 \\
            -@ 12 {input} -o {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"

        rm -rf $target_dir/tmp_{wildcards.BaseName}
        """)