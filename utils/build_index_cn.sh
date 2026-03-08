#!/usr/bin/env bash
# =============================================================================
# build_index_v2.sh - 基因组索引构建脚本 (模板化重构版)
# =============================================================================
# 功能：为多种比对工具构建基因组索引（高代码复用版本）
# 支持工具：bwa-meme, bwa-meth, bwa-mem, astair, batmeth2, bismark-bowtie2,
#           bismark-hisat2, biscuit, bowtie2, bsbolt, bsmapz, fame, gatk,
#           gem3, hisat2, hisat-3n, strobealign, whisper
# 设计目标：在 Docker 容器中运行 (MethFlow)
# 特点：
#   - 模板化函数设计，大幅减少代码重复
#   - 统一的工具配置管理
#   - 支持 -f/--force 强制重新构建
#   - 彩色日志输出和进度显示
# =============================================================================

set -euo pipefail

# =============================================================================
# 颜色定义
# =============================================================================
readonly COLOR_RESET='\033[0m'
readonly COLOR_RED='\033[0;31m'
readonly COLOR_GREEN='\033[0;32m'
readonly COLOR_YELLOW='\033[1;33m'
readonly COLOR_BLUE='\033[0;34m'
readonly COLOR_CYAN='\033[0;36m'
readonly COLOR_MAGENTA='\033[0;35m'

# =============================================================================
# 日志函数
# =============================================================================
log_info() { echo -e "${COLOR_BLUE}[INFO]${COLOR_RESET} $1" >&2; }
log_success() { echo -e "${COLOR_GREEN}[✓]${COLOR_RESET} $1" >&2; }
log_warning() { echo -e "${COLOR_YELLOW}[!]${COLOR_RESET} $1" >&2; }
log_error() { echo -e "${COLOR_RED}[✗]${COLOR_RESET} $1" >&2; }
log_progress() { echo -e "${COLOR_CYAN}[$1/$2]${COLOR_RESET} $3" >&2; }
log_section() {
    echo -e "\n${COLOR_MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${COLOR_RESET}" >&2
    echo -e "${COLOR_MAGENTA}  $1${COLOR_RESET}" >&2
    echo -e "${COLOR_MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${COLOR_RESET}" >&2
}
log_debug() { [[ "${VERBOSE:-false}" == "true" ]] && echo -e "${COLOR_CYAN}[DEBUG]${COLOR_RESET} $1" >&2; }

# =============================================================================
# 变量定义
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METHFLOW_DIR="/opt/MethFlow"
REF=""
OUTPUT_DIR=""
TOOLS="bwa-meme"
THREADS=$(nproc 2>/dev/null || echo 4)
COPY_MODE="symlink"
VERBOSE=false
SKIP_CONDA_CHECK=false
FORCE=false
STROBEALIGN_READ_LENGTH=150

GENOMIC_TOOLS_ENV="/opt/miniforge/envs/genomic_tools"

# =============================================================================
# 工具配置（集中化管理）
# 格式：工具名|conda.yaml路径|二进制路径|构建类型|特殊参数
# =============================================================================
declare -A TOOL_CONFIG=(
    # 需要 conda 环境的工具
    ["bwa-meme"]="conda|rules/bwa-meme/conda.yaml|bwa-meme|single|index -a meme -t {threads} {ref}"
    ["bwa-meth"]="conda|rules/bwa-meth/conda.yaml|bwameth.py|single|index-mem2 {ref}"
    ["bwa-mem"]="conda|rules/astair/conda.yaml|bwa|single|index -a bwtsw -b 500000000 {ref}"
    ["astair"]="conda|rules/astair/conda.yaml|bwa|single|index -a bwtsw -b 500000000 {ref}"
    ["bismark-bowtie2"]="conda|rules/bismark/conda.yaml|bismark_genome_preparation|bismark|--bowtie2 --parallel {half_threads} ."
    ["bismark-hisat2"]="conda|rules/bismark/conda.yaml|bismark_genome_preparation|bismark|--hisat2 --parallel {half_threads} ."
    ["biscuit"]="conda|rules/biscuit/conda.yaml|biscuit|single|index -a bwtsw {ref}"
    ["bowtie2"]="conda|rules/bowtie2/conda.yaml|bowtie2-build|bowtie2|--threads {threads} {ref} {index_prefix}"
    ["bsbolt"]="conda|rules/bsbolt/conda.yaml|bsbolt|bsbolt|Index -G {ref} -DB ."
    ["gatk"]="conda|rules/gatk/conda.yaml|gatk|gatk|CreateSequenceDictionary -R {ref} -O {dict_file}"
    ["strobealign"]="conda|rules/strobealign/conda.yaml|strobealign|strobealign|--create-index -t {threads} -r {read_length} {ref}"
    
    # 不需要 conda 的工具（直接使用容器内版本）
    ["batmeth2"]="system||/opt/MethFlow/resources/BatMeth2/bin/BatMeth2|batmeth2|index -g {ref}"
    ["bsmapz"]="system||special|bsmapz|special"
    ["fame"]="system||/opt/MethFlow/resources/FAME/FAME|fame|--genome {ref} --store_index {ref}.fame"
    ["gem3"]="system||/opt/MethFlow/resources/gem3-mapper/bin/gem-indexer|gem3|-i {ref} -o {index_prefix} -t {threads} -v"
    ["hisat2"]="system||/opt/MethFlow/resources/hisat2/hisat2-build|hisat2|--large-index -p {threads} {ref} {index_prefix}"
    ["hisat-3n"]="system||/opt/MethFlow/resources/hisat-3n/hisat-3n|hisat3n|build"
    ["whisper"]="system||/opt/MethFlow/resources/Whisper/src/whisper-index|whisper|{index_prefix} {ref} {output_dir} {temp_dir}"
)

# 工具索引文件后缀（用于验证）
declare -A TOOL_SUFFIXES=(
    ["bwa-meme"]=".0123 .amb .ann .pac .pos_packed .suffixarray_uint64 .suffixarray_uint64_L0_PARAMETERS .suffixarray_uint64_L1_PARAMETERS .suffixarray_uint64_L2_PARAMETERS"
    ["bwa-meth"]=".bwameth.c2t .bwameth.c2t.0123 .bwameth.c2t.amb .bwameth.c2t.ann .bwameth.c2t.bwt.2bit.64 .bwameth.c2t.pac"
    ["bwa-mem"]=".amb .ann .bwt .pac .sa"
    ["astair"]=".amb .ann .bwt .pac .sa"
    ["batmeth2"]=".batmeth2.fa .batmeth2.fa.amb .batmeth2.fa.ann .batmeth2.fa.bwt .batmeth2.fa.pac .batmeth2.fa.sa .bin .len"
    ["bismark-bowtie2"]="Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"
    ["bismark-hisat2"]="Bisulfite_Genome/CT_conversion/BS_CT.1.ht2 Bisulfite_Genome/GA_conversion/BS_GA.1.ht2"
    ["biscuit"]=".bis.amb .bis.ann .bis.pac .dau.bwt .dau.sa .par.bwt .par.sa"
    ["bowtie2"]=".1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2"
    ["bsbolt"]="BSB_ref.fa BSB_ref.fa.amb BSB_ref.fa.ann BSB_ref.fa.bwt BSB_ref.fa.opac BSB_ref.fa.pac BSB_ref.fa.sa genome_index.pkl"
    ["bsmapz"]=".fa .fai"
    ["fame"]=".fame .fame_strands"
    ["gatk"]=".dict"
    ["gem3"]=".gem .info bs/.gem bs/.info"
    ["hisat2"]=".1.ht2l .2.ht2l .3.ht2l .4.ht2l .5.ht2l .6.ht2l .7.ht2l .8.ht2l"
    ["hisat-3n"]="_c2t.3n.CT.1.ht2 _c2t.3n.GA.1.ht2 _t2c.3n.CT.1.ht2 _t2c.3n.GA.1.ht2"
    ["strobealign"]=".r150.sti"
    ["whisper"]=".whisper_idx.lut_long_dir .whisper_idx.lut_long_rc .whisper_idx.lut_short_dir .whisper_idx.lut_short_rc .whisper_idx.ref_seq_desc .whisper_idx.ref_seq_dir_pck .whisper_idx.ref_seq_rc_pck .whisper_idx.sa_dir .whisper_idx.sa_rc"
)

# =============================================================================
# 帮助信息
# =============================================================================
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

基因组索引构建脚本 (模板化重构版)
设计用于: MethFlow Docker 容器

Required Arguments:
    -r, --ref FILE          参考基因组 FASTA 文件路径
    -o, --output-dir DIR    输出目录

Optional Arguments:
    -t, --tools TOOLS       工具列表，逗号分隔 (默认: bwa-meme)
                            支持: ${!TOOL_CONFIG[*]}
    -m, --methflow-dir DIR  MethFlow 项目根目录
                            (默认: /opt/MethFlow)
    -@, --threads NUM       线程数 (默认: 使用所有 CPU)
    --copy-mode MODE        文件复制方式: copy|symlink (默认: symlink)
    --strobealign-read-length NUM  strobealign 的 reads 长度 (默认: 150)
    -f, --force             强制重新构建索引（即使索引已存在）
    --skip-conda-check      跳过 conda 环境存在性检查
    -v, --verbose           详细输出模式
    -h, --help              显示此帮助信息

Docker 使用示例:
    docker run -it --rm \\
        -v /host/ref:/ref \\
        -v /host/output:/output \\
        MethFlow:latest \\
        /opt/MethFlow/utils/build_index_v2.sh \\
        -r /ref/genome.fa -o /output -t bwa-meme,bwa-meth -@ 16

EOF
}

# =============================================================================
# 基础工具函数
# =============================================================================

# 解析命令行参数
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -r|--ref) REF="$2"; shift 2 ;;
            -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
            -t|--tools) TOOLS="$2"; shift 2 ;;
            -m|--methflow-dir) METHFLOW_DIR="$2"; shift 2 ;;
            -@|--threads) THREADS="$2"; shift 2 ;;
            --copy-mode) COPY_MODE="$2"; shift 2 ;;
            --strobealign-read-length) STROBEALIGN_READ_LENGTH="$2"; shift 2 ;;
            -f|--force) FORCE=true; shift ;;
            --skip-conda-check) SKIP_CONDA_CHECK=true; shift ;;
            -v|--verbose) VERBOSE=true; shift ;;
            -h|--help) usage; exit 0 ;;
            *) log_error "未知参数: $1"; usage; exit 1 ;;
        esac
    done
}

# 验证参数
validate_args() {
    local errors=0
    if [[ -z "$REF" ]]; then
        log_error "缺少必需参数: -r/--ref"
        errors=$((errors + 1))
    elif [[ ! -f "$REF" ]]; then
        log_error "参考基因组文件不存在: $REF"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "缺少必需参数: -o/--output-dir"
        errors=$((errors + 1))
    fi
    
    if [[ "$COPY_MODE" != "copy" && "$COPY_MODE" != "symlink" ]]; then
        log_error "无效的复制模式: $COPY_MODE (必须是 'copy' 或 'symlink')"
        errors=$((errors + 1))
    fi
    
    if [[ $errors -gt 0 ]]; then
        usage
        exit 1
    fi
    
    # 获取绝对路径
    METHFLOW_DIR=$(cd "$METHFLOW_DIR" && pwd)
    mkdir -p "$OUTPUT_DIR"
    OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
    REF=$(realpath "$REF")
}

# 获取 conda 环境路径
get_conda_env() {
    local tool=$1
    local conda_yaml="${TOOL_CONFIG[$tool]#*|}"
    conda_yaml="${conda_yaml%%|*}"
    
    local env_check_script="${SCRIPT_DIR}/get_conda_env_path.sh"
    if [[ ! -f "$env_check_script" ]]; then
        log_error "找不到环境检查脚本: $env_check_script"
        return 1
    fi
    
    local full_conda_yaml="${METHFLOW_DIR}/${conda_yaml}"
    if [[ ! -f "$full_conda_yaml" ]]; then
        log_error "conda.yaml 文件不存在: $full_conda_yaml"
        return 1
    fi
    
    local env_path
    env_path=$("$env_check_script" "$conda_yaml" "$METHFLOW_DIR" 2>/dev/null)
    local exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        local error_msg=$("$env_check_script" "$conda_yaml" "$METHFLOW_DIR" 2>&1 >/dev/null)
        if echo "$error_msg" | grep -q "WARNING"; then
            log_warning "conda 环境尚未创建"
            log_info "提示: 在容器中运行以下命令创建环境:"
            log_info "  cd $METHFLOW_DIR"
            log_info "  snakemake --use-conda --conda-create-envs-only --cores 1"
            if [[ "$SKIP_CONDA_CHECK" == "false" ]]; then
                return 1
            fi
        else
            log_error "无法获取 conda 环境路径"
            return 1
        fi
    fi

    env_path=$(echo "$env_path" | tr -d '\n\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$env_path" ]]; then
        log_error "获取到的 conda 环境路径为空"
        return 1
    fi

    # 检查环境目录是否实际存在
    if [[ ! -d "$env_path" ]]; then
        log_error "conda 环境目录不存在: $env_path"
        log_info "提示: 在容器中运行以下命令创建环境:"
        log_info "  cd $METHFLOW_DIR"
        log_info "  snakemake --use-conda --conda-create-envs-only --cores 1"
        if [[ "$SKIP_CONDA_CHECK" == "false" ]]; then
            return 1
        fi
    fi

    log_debug "conda 环境路径: $env_path"
    printf '%s' "$env_path"
}

# 设置参考基因组（复制或链接）
setup_reference() {
    local ref=$1
    local target_dir=$2
    local ref_name=$(basename "$ref")
    local target_ref="${target_dir}/${ref_name}"
    
    log_info "设置参考基因组: $ref_name"
    log_debug "源文件: $ref"
    log_debug "目标: $target_ref"
    
    # 检查并处理已存在的文件
    if [[ -e "$target_ref" ]]; then
        if [[ "$COPY_MODE" == "symlink" ]] && [[ -L "$target_ref" ]]; then
            local current_link=$(readlink "$target_ref")
            if [[ "$current_link" == "$ref" ]]; then
                log_debug "软链接已存在且指向正确"
                printf '%s' "$target_ref"
                return 0
            else
                log_warning "软链接指向不同源文件，删除旧链接..."
                rm -f "$target_ref" "${target_ref}.fai"
            fi
        elif [[ "$COPY_MODE" == "copy" ]] && [[ -f "$target_ref" ]]; then
            log_warning "目标文件已存在，复用现有文件"
            printf '%s' "$target_ref"
            return 0
        else
            log_warning "目标位置已存在文件，删除旧文件..."
            rm -f "$target_ref" "${target_ref}.fai"
        fi
    fi
    
    # 创建复制或链接
    if [[ "$COPY_MODE" == "copy" ]]; then
        log_info "复制参考基因组到目标目录..."
        cp "$ref" "$target_ref"
        [[ -f "${ref}.fai" ]] && cp "${ref}.fai" "${target_ref}.fai"
    else
        log_info "创建参考基因组软链接..."
        ln -s "$ref" "$target_ref"
        [[ -f "${ref}.fai" ]] && ln -s "${ref}.fai" "${target_ref}.fai"
    fi
    
    log_success "参考基因组已设置"
    printf '%s' "$target_ref"
}

# 检查索引是否存在
check_index_exists() {
    local tool=$1
    local ref=$2
    local suffixes="${TOOL_SUFFIXES[$tool]}"
    local output_dir=$(dirname "$ref")
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    
    local missing=0
    for suffix in $suffixes; do
        local check_file
        if [[ "$suffix" == */* ]]; then
            # 复合路径（如 Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 或 bs/.gem）
            if [[ "$tool" == "gem3" ]]; then
                # gem3 的 bs 子目录中的文件：bs/.gem → bs/GRCh38.gem
                local dir_part="${suffix%/*}"   # bs
                local file_part="${suffix##*/}" # .gem
                check_file="${output_dir}/${dir_part}/${index_prefix}${file_part}"
            else
                check_file="${output_dir}/${suffix}"
            fi
        elif [[ "$tool" == "whisper" || "$tool" == "gatk" || "$tool" == "bowtie2" || "$tool" == "hisat2" || "$tool" == "hisat-3n" ]]; then
            # 这些工具的索引文件使用 ${index_prefix}${suffix} 格式（如 GRCh38.dict, GRCh38.1.bt2, GRCh38_c2t.3n.CT.1.ht2）
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "gem3" ]]; then
            # gem3 普通索引文件使用 ${index_prefix}${suffix} 格式（如 GRCh38.gem）
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "bsbolt" ]]; then
            # bsbolt 索引文件直接使用 suffix 作为文件名（如 BSB_ref.fa, genome_index.pkl）
            check_file="${output_dir}/${suffix}"
        elif [[ "$tool" == "bsmapz" ]]; then
            # bsmapz 的 .fa 后缀就是参考文件本身，不是附加后缀
            if [[ "$suffix" == ".fa" ]]; then
                check_file="${output_dir}/${ref_name}"
            else
                check_file="${output_dir}/${ref_name}${suffix}"
            fi
        else
            check_file="${ref}${suffix}"
        fi

        [[ ! -f "$check_file" ]] && ((missing++))
    done

    return $missing
}

# =============================================================================
# 通用模板函数
# =============================================================================

# 通用索引构建函数
build_index_generic() {
    local tool=$1
    local ref=$2
    local output_dir=$3
    local threads=$4
    local conda_env=$5
    
    log_info "开始构建 $tool 索引..."
    log_info "参考基因组: $(basename "$ref")"
    log_info "线程数: $threads"
    log_debug "工作目录: $output_dir"
    
    cd "$output_dir"
    
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    local half_threads=$((threads / 2))
    [[ $half_threads -lt 1 ]] && half_threads=1
    
    # 获取工具配置
    # 格式: env_type|conda_yaml|binary|build_type|args
    local config="${TOOL_CONFIG[$tool]}"
    local env_type="${config%%|*}"
    config="${config#*|}"
    local conda_yaml="${config%%|*}"
    config="${config#*|}"
    local binary="${config%%|*}"
    config="${config#*|}"
    local build_type="${config%%|*}"
    local args="${config#*|}"
    
    # 特殊处理：hisat-3n 需要构建两次
    if [[ "$tool" == "hisat-3n" ]]; then
        build_hisat3n_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # 特殊处理：gem3 需要构建两次（普通 + bisulfite）
    if [[ "$tool" == "gem3" ]]; then
        build_gem3_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # 特殊处理：batmeth2 需要构建两次（默认 + RRBS）
    if [[ "$tool" == "batmeth2" ]]; then
        build_batmeth2_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # 特殊处理：bsmapz 不需要构建索引
    if [[ "$tool" == "bsmapz" ]]; then
        build_bsmapz_ref "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # 特殊处理：whisper
    if [[ "$tool" == "whisper" ]]; then
        build_whisper_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # 替换变量到 args
    args="${args//\{threads\}/$threads}"
    args="${args//\{half_threads\}/$half_threads}"
    args="${args//\{ref\}/$ref}"
    args="${args//\{index_prefix\}/$index_prefix}"
    args="${args//\{read_length\}/$STROBEALIGN_READ_LENGTH}"
    args="${args//\{dict_file\}/${index_prefix}.dict}"
    
    log_info "运行: $binary $args"

    # 执行命令
    local full_cmd
    if [[ -n "$conda_env" ]]; then
        full_cmd="mamba run -p '$conda_env' $binary $args"
    else
        full_cmd="$binary $args"
    fi
    
    if [[ "$VERBOSE" == "true" ]]; then
        eval "$full_cmd" 2>&1 | while IFS= read -r line; do
            echo -e "${COLOR_CYAN}    [$tool]${COLOR_RESET} $line"
        done
    else
        eval "$full_cmd" 2>&1 | while IFS= read -r line; do
            if echo "$line" | grep -qiE "(index|build|done|finish|complete|step|process|creating)"; then
                echo -e "${COLOR_CYAN}    [$tool]${COLOR_RESET} $line"
            fi
        done
    fi
    
    local exit_code=${PIPESTATUS[0]}
    
    if [[ $exit_code -eq 0 ]]; then
        log_success "$tool 索引构建完成"
        return 0
    else
        log_error "$tool 索引构建失败 (退出码: $exit_code)"
        return 1
    fi
}

# 通用索引验证函数
verify_index_generic() {
    local tool=$1
    local ref=$2
    local output_dir=$3
    
    log_info "验证 $tool 索引文件..."
    
    local suffixes="${TOOL_SUFFIXES[$tool]}"
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    
    local count=0
    local total=0
    
    for suffix in $suffixes; do
        ((total++))
        local check_file
        if [[ "$suffix" == */* ]]; then
            # 复合路径（如 Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 或 bs/.gem）
            if [[ "$tool" == "gem3" ]]; then
                # gem3 的 bs 子目录中的文件：bs/.gem → bs/GRCh38.gem
                local dir_part="${suffix%/*}"   # bs
                local file_part="${suffix##*/}" # .gem
                check_file="${output_dir}/${dir_part}/${index_prefix}${file_part}"
            else
                check_file="${output_dir}/${suffix}"
            fi
        elif [[ "$tool" == "whisper" || "$tool" == "gatk" || "$tool" == "bowtie2" || "$tool" == "hisat2" || "$tool" == "hisat-3n" ]]; then
            # 这些工具的索引文件使用 ${index_prefix}${suffix} 格式（如 GRCh38.dict, GRCh38.1.bt2, GRCh38_c2t.3n.CT.1.ht2）
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "gem3" ]]; then
            # gem3 普通索引文件使用 ${index_prefix}${suffix} 格式（如 GRCh38.gem）
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "bsbolt" ]]; then
            # bsbolt 索引文件直接使用 suffix 作为文件名（如 BSB_ref.fa, genome_index.pkl）
            check_file="${output_dir}/${suffix}"
        elif [[ "$tool" == "bsmapz" ]]; then
            # bsmapz 的 .fa 后缀就是参考文件本身，不是附加后缀
            if [[ "$suffix" == ".fa" ]]; then
                check_file="${output_dir}/${ref_name}"
            else
                check_file="${output_dir}/${ref_name}${suffix}"
            fi
        else
            check_file="${ref}${suffix}"
        fi

        if [[ -f "$check_file" ]]; then
            local size=$(du -h "$check_file" 2>/dev/null | cut -f1)
            log_success "  $(basename "$check_file") (${size})"
            ((count++))
        else
            log_warning "  $(basename "$check_file") 缺失"
        fi
    done

    if [[ $count -eq $total ]]; then
        log_success "索引验证通过 ($count/$total 文件存在)"
        return 0
    elif [[ $count -ge $((total - 2)) ]]; then
        log_warning "索引部分存在 ($count/$total 文件存在)，可能可用"
        return 0
    else
        log_error "索引验证失败 ($count/$total 文件存在)"
        return 1
    fi
}

# =============================================================================
# 特殊工具处理函数（需要多步构建的）
# =============================================================================

build_hisat3n_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local hisat3n_bin="/opt/MethFlow/resources/hisat-3n/hisat-3n"

    log_info "开始构建 hisat-3n 索引..."
    log_info "参考基因组: $(basename "$ref")"
    log_info "线程数: $threads"
    log_debug "工作目录: $output_dir"
    log_debug "hisat-3n 路径: $hisat3n_bin"

    # 检查 hisat-3n 是否存在
    if [[ ! -f "$hisat3n_bin" ]]; then
        log_error "hisat-3n 二进制文件不存在: $hisat3n_bin"
        return 1
    fi

    # 切换到工作目录
    cd "$output_dir"

    local ref_name
    ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"

    # 检查索引是否已存在
    if [[ "$FORCE" == "false" ]] && check_index_exists "hisat-3n" "$ref"; then
        log_warning "hisat-3n 索引已存在，跳过构建（使用 -f/--force 强制重新构建）"
        return 0
    fi

    # 步骤 1: 构建 C2T 索引 (--base-change C,T)
    log_info "步骤 1/2: 构建 C2T 索引..."
    log_info "运行: hisat-3n-build --base-change C,T -p ${threads} $(basename "$ref") ${index_prefix}_c2t"

    if [[ "$VERBOSE" == "true" ]]; then
        "$hisat3n_bin"-build --base-change C,T -p "$threads" "$ref" "${index_prefix}_c2t" 2>&1 | while IFS= read -r line; do
            echo -e "${COLOR_CYAN}    [hisat-3n-c2t]${COLOR_RESET} $line"
        done
    else
        "$hisat3n_bin"-build --base-change C,T -p "$threads" "$ref" "${index_prefix}_c2t" 2>&1 | while IFS= read -r line; do
            if echo "$line" | grep -qiE "(index|build|done|finish|complete|step|process|settings)"; then
                echo -e "${COLOR_CYAN}    [hisat-3n-c2t]${COLOR_RESET} $line"
            fi
        done
    fi

    local exit_code=${PIPESTATUS[0]}
    if [[ $exit_code -ne 0 ]]; then
        log_error "hisat-3n C2T 索引构建失败 (退出码: $exit_code)"
        return 1
    fi
    log_success "C2T 索引构建完成"

    # 步骤 2: 构建 T2C 重复索引 (--base-change T,C --repeat-index)
    log_info "步骤 2/2: 构建 T2C 重复索引..."
    log_info "运行: hisat-3n-build --base-change T,C --repeat-index -p ${threads} $(basename "$ref") ${index_prefix}_t2c"

    if [[ "$VERBOSE" == "true" ]]; then
        "$hisat3n_bin"-build --base-change T,C --repeat-index -p "$threads" "$ref" "${index_prefix}_t2c" 2>&1 | while IFS= read -r line; do
            echo -e "${COLOR_CYAN}    [hisat-3n-t2c]${COLOR_RESET} $line"
        done
    else
        "$hisat3n_bin"-build --base-change T,C --repeat-index -p "$threads" "$ref" "${index_prefix}_t2c" 2>&1 | while IFS= read -r line; do
            if echo "$line" | grep -qiE "(index|build|done|finish|complete|step|process|settings)"; then
                echo -e "${COLOR_CYAN}    [hisat-3n-t2c]${COLOR_RESET} $line"
            fi
        done
    fi

    exit_code=${PIPESTATUS[0]}
    if [[ $exit_code -ne 0 ]]; then
        log_error "hisat-3n T2C 重复索引构建失败 (退出码: $exit_code)"
        return 1
    fi
    log_success "T2C 重复索引构建完成"

    log_success "hisat-3n 所有索引构建完成"
    return 0
}

build_gem3_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local gem3_bin="/opt/MethFlow/resources/gem3-mapper/bin/gem-indexer"
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    
    log_info "步骤 1/2: 构建普通索引..."
    "$gem3_bin" -i "$ref" -o "${output_dir}/${index_prefix}" -t "$threads" -v 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [gem3]${COLOR_RESET} $line"; done
    
    log_info "步骤 2/2: 构建 bisulfite 索引..."
    local bs_dir="${output_dir}/bs"
    mkdir -p "$bs_dir"
    ln -sf "$ref" "${bs_dir}/${ref_name}" 2>/dev/null || true
    cd "$bs_dir"
    "$gem3_bin" -i "$ref" -o "${index_prefix}" -b -t "$threads" -v 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [gem3-bs]${COLOR_RESET} $line"; done
    cd "$output_dir"
    
    log_success "gem3 索引构建完成"
}

build_batmeth2_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local batmeth2_bin="/opt/MethFlow/resources/BatMeth2/bin/BatMeth2"
    
    log_info "步骤 1/2: 构建默认索引..."
    "$batmeth2_bin" index -g "$ref" 2>&1 | grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [BatMeth2]${COLOR_RESET} $line"; done
    
    log_info "步骤 2/2: 构建 RRBS 索引..."
    local rrbs_dir="${output_dir}/rrbs"
    mkdir -p "$rrbs_dir"
    ln -sf "$ref" "${rrbs_dir}/$(basename "$ref")" 2>/dev/null || true
    cd "$rrbs_dir"
    "$batmeth2_bin" index_rrbs -g "$(basename "$ref")" 2>&1 | grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [BatMeth2-RRBS]${COLOR_RESET} $line"; done
    cd "$output_dir"
    
    log_success "BatMeth2 索引构建完成"
}

build_bsmapz_ref() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    
    log_info "准备 bsmapz 参考基因组..."
    local ref_name=$(basename "$ref")
    local target_ref="${output_dir}/${ref_name}"
    
    # 检查行长度
    local line_len=$(grep -v '^>' "$ref" | grep -v '^$' | head -1 | wc -c)
    if [[ $line_len -le 70 ]]; then
        log_info "行长度符合要求，创建软链接..."
        ln -sf "$ref" "$target_ref"
    else
        log_info "行长度超过 70，重新格式化..."
        mamba run -p "$GENOMIC_TOOLS_ENV" seqtk seq -l 70 "$ref" > "$target_ref"
    fi
    
    # 生成 .fai
    if [[ ! -f "${target_ref}.fai" ]]; then
        mamba run -p "$GENOMIC_TOOLS_ENV" samtools faidx "$target_ref"
    fi
    
    log_success "bsmapz 参考基因组准备完成"
}

build_whisper_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local whisper_bin="/opt/MethFlow/resources/Whisper/src/whisper-index"
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    local temp_dir="${output_dir}/temp"
    
    mkdir -p "$temp_dir"
    log_info "运行: whisper-index $index_prefix $ref_name $output_dir $temp_dir"
    
    "$whisper_bin" "$index_prefix" "$ref" "$output_dir" "$temp_dir" 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [whisper]${COLOR_RESET} $line"; done
    
    rm -rf "$temp_dir"
    log_success "whisper 索引构建完成"
}

# =============================================================================
# 主流程函数
# =============================================================================

# 处理单个工具
process_tool() {
    local tool=$1
    local current=$2
    local total=$3
    
    log_section "处理工具: $tool"
    log_progress "$current" "$total" "开始处理 $tool"
    
    # 检查工具是否支持
    if [[ -z "${TOOL_CONFIG[$tool]:-}" ]]; then
        log_error "不支持的工具: $tool"
        return 1
    fi
    
    # 解析工具配置
    local config="${TOOL_CONFIG[$tool]}"
    local env_type="${config%%|*}"
    config="${config#*|}"
    
    # 获取 conda 环境（如果需要）
    local conda_env=""
    if [[ "$env_type" == "conda" ]]; then
        conda_env=$(get_conda_env "$tool")
        if [[ $? -ne 0 ]]; then
            log_error "无法获取 $tool 的 conda 环境"
            return 1
        fi
        log_success "找到 conda 环境"
    else
        log_info "工具 $tool 使用容器内安装的版本，跳过 conda 环境"
    fi
    
    # 创建工具输出目录
    local tool_dir="${OUTPUT_DIR}/${tool}"
    mkdir -p "$tool_dir"
    log_info "工具输出目录: $tool_dir"
    
    # 设置参考基因组
    local target_ref
    target_ref=$(setup_reference "$REF" "$tool_dir")
    
    # 检查索引是否已存在
    if [[ "$FORCE" == "false" ]] && check_index_exists "$tool" "$target_ref"; then
        log_warning "$tool 索引已存在，跳过构建（使用 -f/--force 强制重新构建）"
        verify_index_generic "$tool" "$target_ref" "$tool_dir"
        return 0
    fi
    
    # 构建索引
    local build_success=false
    build_index_generic "$tool" "$target_ref" "$tool_dir" "$THREADS" "$conda_env"
    
    if [[ $? -eq 0 ]]; then
        verify_index_generic "$tool" "$target_ref" "$tool_dir"
        if [[ $? -eq 0 ]]; then
            build_success=true
        fi
    fi
    
    if [[ "$build_success" == "true" ]]; then
        log_success "$tool 索引构建成功"
        return 0
    else
        log_error "$tool 索引构建失败"
        return 1
    fi
}

# 显示系统信息
show_system_info() {
    log_section "系统信息"
    if [[ -f /.dockerenv ]] || grep -q docker /proc/1/cgroup 2>/dev/null; then
        log_info "运行环境: Docker 容器"
    else
        log_info "运行环境: 宿主机"
    fi
    log_info "CPU 核心数: $(nproc)"
    log_info "可用内存: $(free -h 2>/dev/null | awk '/^Mem:/ {print $2}' || echo "N/A")"
    log_info "当前用户: $(whoami)"
    log_info "工作目录: $(pwd)"
}

# =============================================================================
# 主函数
# =============================================================================
main() {
    log_section "基因组索引构建工具 (模板化重构版)"
    
    # 解析参数
    parse_args "$@"
    validate_args
    
    # 显示系统信息
    show_system_info
    
    # 显示配置信息
    log_section "配置信息"
    log_info "参考基因组: $REF"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "工具列表: $TOOLS"
    log_info "Methflow 目录: $METHFLOW_DIR"
    log_info "线程数: $THREADS"
    log_info "复制模式: $COPY_MODE"
    log_info "强制重新构建: $FORCE"
    log_info "strobealign reads 长度: ${STROBEALIGN_READ_LENGTH}bp"
    log_info "详细模式: $VERBOSE"
    
    # 分割工具列表
    IFS=',' read -ra TOOL_ARRAY <<< "$TOOLS"
    local total=${#TOOL_ARRAY[@]}
    log_info "共需处理 $total 个工具"
    
    # 处理每个工具
    local success_count=0
    local fail_count=0
    local i=1
    
    for tool in "${TOOL_ARRAY[@]}"; do
        tool=$(echo "$tool" | xargs)
        [[ -z "$tool" ]] && continue
        
        if process_tool "$tool" "$i" "$total"; then
            ((success_count++))
        else
            ((fail_count++))
        fi
        
        ((i++))
        echo ""
    done
    
    # 显示汇总
    log_section "构建完成"
    log_info "总计: $total 个工具"
    log_success "成功: $success_count 个"
    
    if [[ $fail_count -gt 0 ]]; then
        log_error "失败: $fail_count 个"
        exit 1
    fi
    
    log_success "所有索引构建完成！"
    log_info "输出目录: $OUTPUT_DIR"
    exit 0
}

# 执行主函数
main "$@"

