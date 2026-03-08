configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
pufferfish index -r ./BL.fa -o . -p 16
pufferfish align \
    -i /hot_warm_data/ref/quartet/DNA/custom_genome/BL/pufferfish \
    --mate1 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    --mate2 /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz \
    -o /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/pufferfish.sam \
    --threads 16 \
    --allowSoftclip --allowOverhangSoftclip \
    --writeQualities --verbose --genomicReads
"""


rule pufferfish_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/pufferfish/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/pufferfish/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/pufferfish/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["pufferfish"][wildcards.BaseName.split('_')[1]],
        extra_params = config["pufferfish"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        pufferfish_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/pufferfish"

        resources/pufferfish/build/src/pufferfish align \\
            -i {params.ref} --mate1 {input.r1} --mate2 {input.r2} \\
            -o ${{pufferfish_dir}}/pufferfish.sam \\
            --threads {threads} {params.extra_params}
        samtools sort \\
            -@ {threads} \\
            ${{pufferfish_dir}}/pufferfish.sam |\\
        samtools addreplacerg - \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} \\
            -o {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
    
        rm ${{pufferfish_dir}}/pufferfish.sam
        """)
