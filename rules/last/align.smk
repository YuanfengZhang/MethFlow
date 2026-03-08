configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
export PATH=/mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/src:$PATH
export PATH=/mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/script:$PATH

lastdb -c -u BISF -P 64 BL_BISF BL.fa
lastdb -c -u BISR -P 64 BL_BISR BL.fa &
lastdb -c -u NEAR -P 64 BL_NEAR BL.fa & wait

time fastq-interleave \
   /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.R1.fq.gz \
   /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.R2.fq.gz |\
   > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq

time seqtk sample \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq \
    20000 \
    > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM_n20000.fq

# no -P param nor mbuffer because of nearly no improvement
time lastal \
    -Q1 -D 1000 -i 1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISF \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM_n20000.fq |\
    last-pair-probs -e > \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISF.probs
BISF_M=$(grep \
            "estimated mean distance" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISF.probs |\
            awk -F': ' '{print $2}')
BISF_S=$(grep \
            "estimated standard deviation" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISF.probs |\
            awk -F': ' '{print $2}')

time lastal \
    -Q1 -D 1000 -i 1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISR \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM_n20000.fq |\
    last-pair-probs -e > \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISR.probs
BISR_M=$(grep \
            "estimated mean distance" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISR.probs |\
            awk -F': ' '{print $2}')
BISR_S=$(grep \
            "estimated standard deviation" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISR.probs |\
            awk -F': ' '{print $2}')

time lastal \
    -Q1 -D 1000 -i 1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_NEAR \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM_n20000.fq |\
    last-pair-probs -e > \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.NEAR.probs
NEAR_M=$(grep \
            "estimated mean distance" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.NEAR.probs |\
            awk -F': ' '{print $2}')
NEAR_S=$(grep \
            "estimated standard deviation" \
            /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.NEAR.probs |\
            awk -F': ' '{print $2}')

last-train \
    --sample-number 20000 -P 64 -X 1 -Q1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISF \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq \
    > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_BISF.train
last-train \
    --sample-number 20000 -P 64 -X 1 -Q1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISR \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq \
    > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_BISR.train &

last-train \
    --sample-number 20000 -P 64 -X 1 -Q1 \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_NEAR \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq \
    > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_NEAR.train & wait

parallel \
    --gnu \
    --citation \
    --pipe \
    -L64 -j64 \
    "lastal \
    -Q1 -i1 -C 2 \
    -p /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_BISF.train \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISF |\
     last-split -m 0.9 -d 2 -n -fMAF+ |\
     /mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/src/last-split-pe \
     -f ${BISF_M} -s ${BISF_S} " \
     < /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq |\
     mbuffer -m 4G -q |\
     samtools sort -n - |\
     > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISF.sam

parallel \
    --gnu \
    --pipe \
    -L64 -j64 \
    "lastal \
    -Q1 -i1 -C 2 \
    -p /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_BISR.train \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_BISR |\
     last-split -m 0.9 -d 2 -n -fMAF+ |\
     /mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/src/last-split-pe \
     -f ${BISR_M} -s ${BISR_S} " \
     < /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq |\
     mbuffer -m 4G -q \
     > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISR.sam

parallel \
    --gnu \
    --pipe \
    -L64 -j64 \
    "lastal \
    -Q1 -i1 -C 2 \
    -p /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM_NEAR.train \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL_NEAR |\
     last-split -m 0.9 -d 2 -n -fMAF+ |\
     /mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/src/last-split-pe \
     -f ${NEAR_M} -s ${NEAR_S} " \
     < /mnt/eqa/zhangyuanfeng/methylation/MethFlow/input/EM_HF_1_IPM.fq |\
     mbuffer -m 4G -q \
     > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.NEAR.sam

python \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/resources/last-split-pe/scripts/printSamHeader.py \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/last/BL.fa \
    > /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/header.sam

python \
    /mnt/eqa/zhangyuanfeng/methylation/MethFlow/rules/last/merge.py \
    --bisf /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISF.sam \
    --bisr /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.BISR.sam \
    --near /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/EM_HF_1_IPM.NEAR.sam \
    --sam_tmp /mnt/eqa/zhangyuanfeng/methylation/MethFlow/test/align/last/sam.tmp

cat header.sam sam.tmp | samtools sort --output-fmt bam,level=9 > EM_HF_1_IPM.bam

"""


rule last_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/last/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/last/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/last/{BaseName}.align.benchmark"
    params:
        bisf_db      = lambda wildcards: config["ref"]["last"][wildcards.BaseName.split('_')[1]]["BISF"],
        bisr_db      = lambda wildcards: config["ref"]["last"][wildcards.BaseName.split('_')[1]]["BISR"],
        near_db      = lambda wildcards: config["ref"]["last"][wildcards.BaseName.split('_')[1]]["NEAR"],
        ref          = lambda wildcards: config["ref"]["last"][wildcards.BaseName.split('_')[1]]["fasta"],
        extra_params = config["last"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        last_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/last"

        mkdir -p $last_dir

        fastq-interleave \\
            {input.r1} {input.r2} \\
            > ${{last_dir}}/{wildcards.BaseName}.fq

        seqtk sample \\
            ${{last_dir}}/{wildcards.BaseName}.fq \\
            20000 \\
            > ${{last_dir}}/{wildcards.BaseName}_20000.fq
        
        case "$LIB" in
            BS|EM|RR)
                lastal \\
                    -Q1 -D 1000 -i 1 \\
                    {params.bisf_db} \\
                    ${{last_dir}}/{wildcards.BaseName}_20000.fq |\\
                last-pair-probs -e > \\
                    ${{last_dir}}/{wildcards.BaseName}.BISF.probs
                BISF_M=$(grep \\
                            "estimated mean distance" \
                            ${{last_dir}}/{wildcards.BaseName}.BISF.probs |\\
                            awk -F': ' '{{print $2}}')
                BISF_S=$(grep \\
                            "estimated standard deviation" \\
                            ${{last_dir}}/{wildcards.BaseName}.BISF.probs |\\
                            awk -F': ' '{{print $2}}')
                lastal \\
                    -Q1 -D 1000 -i 1 \\
                    {params.bisr_db} \\
                    ${{last_dir}}/{wildcards.BaseName}_20000.fq |\\
                last-pair-probs -e > \\
                    ${{last_dir}}/{wildcards.BaseName}.BISR.probs
                BISR_M=$(grep \\
                            "estimated mean distance" \\
                            ${{last_dir}}/{wildcards.BaseName}.BISR.probs |\\
                            awk -F': ' '{{print $2}}')
                BISR_S=$(grep \\
                            "estimated standard deviation" \\
                            ${{last_dir}}/{wildcards.BaseName}.BISR.probs |\\
                            awk -F': ' '{{print $2}}')
                ;;
            *)
                ;;
        esac

        lastal \\
            -Q1 -D 1000 -i 1 \\
            {params.near_db} \\
            ${{last_dir}}/{wildcards.BaseName}_20000.fq |\\
        last-pair-probs -e > \\
            ${{last_dir}}/{wildcards.BaseName}.NEAR.probs
        NEAR_M=$(grep \\
                    "estimated mean distance" \\
                    ${{last_dir}}/{wildcards.BaseName}.NEAR.probs |\\
                    awk -F': ' '{{print $2}}')
        NEAR_S=$(grep \\
                    "estimated standard deviation" \\
                    ${{last_dir}}/{wildcards.BaseName}.NEAR.probs |\\
                    awk -F': ' '{{print $2}}')

        case "$LIB" in
            BS|EM|RR)
                last-train \\
                    --sample-number 20000 -P {threads} -X 1 -Q1 \\
                    {params.bisf_db} ${{last_dir}}/{wildcards.BaseName}.fq \\
                    > ${{last_dir}}/{wildcards.BaseName}.BISF.train
                last-train \\
                    --sample-number 20000 -P {threads} -X 1 -Q1 \\
                    {params.bisr_db} ${{last_dir}}/{wildcards.BaseName}.fq \\
                    > ${{last_dir}}/{wildcards.BaseName}.BISR.train
                ;;
            *)
                ;;
        esac

        last-train \\
            --sample-number 20000 -P {threads} -X 1 -Q1 \\
            {params.near_db} ${{last_dir}}/{wildcards.BaseName}.fq \\
            > ${{last_dir}}/{wildcards.BaseName}.NEAR.train

        case "$LIB" in
            BS|EM|RR)
                parallel \\
                    --gnu \\
                    --pipe \\
                    -L64 -j{threads} \\
                    "lastal \\
                    -Q1 -i1 -C 2 \\
                    -p ${{last_dir}}/{wildcards.BaseName}.BISF.train {params.bisf_db} |\\
                    last-split -m 0.9 -d 2 -n -fMAF+ |\\
                    resources/last-split-pe/src/last-split-pe \\
                    -f ${{BISF_M}} -s ${{BISF_S}} " \\
                    < ${{last_dir}}/{wildcards.BaseName}.fq \\
                    > ${{last_dir}}/{wildcards.BaseName}.BISF.sam

                parallel \\
                    --gnu \\
                    --pipe \\
                    -L64 -j{threads} \\
                    "lastal \\
                    -Q1 -i1 -C 2 \\
                    -p ${{last_dir}}/{wildcards.BaseName}.BISR.train {params.bisr_db} |\\
                    last-split -m 0.9 -d 2 -n -fMAF+ |\\
                    resources/last-split-pe/src/last-split-pe \\
                    -f ${{BISR_M}} -s ${{BISR_S}} " \\
                    < ${{last_dir}}/{wildcards.BaseName}.fq \\
                    > ${{last_dir}}/{wildcards.BaseName}.BISR.sam
                ;;
            *)
                ;;
        esac

        parallel \\
            --gnu \\
            --pipe \\
            -L64 -j{threads} \\
            "lastal \\
            -Q1 -i1 -C 2 \\
            -p ${{last_dir}}/{wildcards.BaseName}.NEAR.train {params.near_db} |\\
            last-split -m 0.9 -d 2 -n -fMAF+ |\\
            resources/last-split-pe/src/last-split-pe \\
            -f ${{NEAR_M}} -s ${{NEAR_S}} " \\
            < ${{last_dir}}/{wildcards.BaseName}.fq \\
            > ${{last_dir}}/{wildcards.BaseName}.NEAR.sam
        
        python \\
            resources/last-split-pe/scripts/printSamHeader.py \\
            {params.ref} \\
            > ${{last_dir}}/header.sam
        case "$LIB" in
            BS|EM|RR)
                python \\
                    rules/last/merge.py \
                    --bisf ${{last_dir}}/{wildcards.BaseName}.BISF.sam \\
                    --bisr ${{last_dir}}/{wildcards.BaseName}.BISR.sam \\
                    --near ${{last_dir}}/{wildcards.BaseName}.NEAR.sam \\\
                    --sam_tmp ${{last_dir}}/sam.tmp
                ;;
            *)
                mv ${{last_dir}}/{wildcards.BaseName}.NEAR.sam \\
                    ${{last_dir}}/sam.tmp
                ;;
        esac

        cat ${{last_dir}}/header.sam ${{last_dir}}/sam.tmp |\\
        samtools sort -@ {threads} - |\\
        samtools addreplacerg - \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"

        cd ${{last_dir}}
        rm *.fq *.sam *.probs *.train *.tmp
        """)
