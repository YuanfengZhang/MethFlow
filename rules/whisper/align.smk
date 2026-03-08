configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
whisper-index BL BL.fa . ./temp
whisper -rp \
    -t 24 -store-BAM -gzipped-SAM 9 \
    -out /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/whisper \
    -temp ./temp \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/Whisper/BL \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R1.fq.gz \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/BS_HF_1_MX.R2.fq.gz
"""


rule whisper_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/whisper/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/whisper/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/whisper/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["whisper"][wildcards.BaseName.split('_')[1]],
        extra_params = config["whisper"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        whisper_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/whisper"

        mkdir "${{whisper_dir}}/whisper_tmp"

        resources/Whisper/src/whisper -rp \\
            -t {threads} {params.extra_params} \\
            -out ${{whisper_dir}}/whisper \\
            -temp "${{whisper_dir}}/whisper_tmp" \\
            {params.ref} {input.r1} {input.r2}

        samtools sort \\
            --output-fmt bam,level=9 \\
            -@ {threads} \\
            ${{whisper_dir}}/whisper.bam |\\
        samtools addreplacerg - \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} \\
            -o {output.bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
    
        rm -rf ${{whisper_dir}}/whisper.bam ${{whisper_dir}}/whisper_tmp
        """)
