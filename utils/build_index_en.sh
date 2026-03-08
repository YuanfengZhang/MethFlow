#!/usr/bin/env bash
# =============================================================================
# build_index_en.sh - Genome Index Building Script (Template-based Refactored Version)
# =============================================================================
# Function: Build genome indices for multiple alignment tools (high code reuse version)
# Supported tools: bwa-meme, bwa-meth, bwa-mem, astair, batmeth2, bismark-bowtie2,
#                  bismark-hisat2, biscuit, bowtie2, bsbolt, bsmapz, fame, gatk,
#                  gem3, hisat2, hisat-3n, strobealign, whisper
# Design goal: Run inside Docker container (MethFlow)
# Features:
#   - Template-based function design, greatly reducing code duplication
#   - Unified tool configuration management
#   - Support for -f/--force to force rebuild
#   - Colored log output and progress display
# =============================================================================

set -euo pipefail

# =============================================================================
# Color Definitions
# =============================================================================
readonly COLOR_RESET='\033[0m'
readonly COLOR_RED='\033[0;31m'
readonly COLOR_GREEN='\033[0;32m'
readonly COLOR_YELLOW='\033[1;33m'
readonly COLOR_BLUE='\033[0;34m'
readonly COLOR_CYAN='\033[0;36m'
readonly COLOR_MAGENTA='\033[0;35m'

# =============================================================================
# Logging Functions
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
# Variable Definitions
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
# Tool Configuration (Centralized Management)
# Format: tool_name|conda.yaml_path|binary_path|build_type|special_args
# =============================================================================
declare -A TOOL_CONFIG=(
    # Tools requiring conda environment
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
    
    # Tools not requiring conda (using container-installed versions directly)
    ["batmeth2"]="system||/opt/MethFlow/resources/BatMeth2/bin/BatMeth2|batmeth2|index -g {ref}"
    ["bsmapz"]="system||special|bsmapz|special"
    ["fame"]="system||/opt/MethFlow/resources/FAME/FAME|fame|--genome {ref} --store_index {ref}.fame"
    ["gem3"]="system||/opt/MethFlow/resources/gem3-mapper/bin/gem-indexer|gem3|-i {ref} -o {index_prefix} -t {threads} -v"
    ["hisat2"]="system||/opt/MethFlow/resources/hisat2/hisat2-build|hisat2|--large-index -p {threads} {ref} {index_prefix}"
    ["hisat-3n"]="system||/opt/MethFlow/resources/hisat-3n/hisat-3n|hisat3n|build"
    ["whisper"]="system||/opt/MethFlow/resources/Whisper/src/whisper-index|whisper|{index_prefix} {ref} {output_dir} {temp_dir}"
)

# Tool index file suffixes (for verification)
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
# Help Information
# =============================================================================
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Genome Index Building Script (Template-based Refactored Version)
Designed for: MethFlow Docker container

Required Arguments:
    -r, --ref FILE          Reference genome FASTA file path
    -o, --output-dir DIR    Output directory

Optional Arguments:
    -t, --tools TOOLS       Tool list, comma-separated (default: bwa-meme)
                            Supported: ${!TOOL_CONFIG[*]}
    -m, --methflow-dir DIR  MethFlow project root directory
                            (default: /opt/MethFlow)
    -@, --threads NUM       Number of threads (default: use all CPUs)
    --copy-mode MODE        File copy method: copy|symlink (default: symlink)
    --strobealign-read-length NUM  Read length for strobealign (default: 150)
    -f, --force             Force rebuild indices (even if they already exist)
    --skip-conda-check      Skip conda environment existence check
    -v, --verbose           Verbose output mode
    -h, --help              Show this help message

Docker Usage Example:
    docker run -it --rm \\
        -v /host/ref:/ref \\
        -v /host/output:/output \\
        MethFlow:latest \\
        /opt/MethFlow/utils/build_index_en.sh \\
        -r /ref/genome.fa -o /output -t bwa-meme,bwa-meth -@ 16

EOF
}

# =============================================================================
# Basic Utility Functions
# =============================================================================

# Parse command line arguments
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
            *) log_error "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done
}

# Validate arguments
validate_args() {
    local errors=0
    if [[ -z "$REF" ]]; then
        log_error "Missing required argument: -r/--ref"
        errors=$((errors + 1))
    elif [[ ! -f "$REF" ]]; then
        log_error "Reference genome file does not exist: $REF"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "Missing required argument: -o/--output-dir"
        errors=$((errors + 1))
    fi
    
    if [[ "$COPY_MODE" != "copy" && "$COPY_MODE" != "symlink" ]]; then
        log_error "Invalid copy mode: $COPY_MODE (must be 'copy' or 'symlink')"
        errors=$((errors + 1))
    fi
    
    if [[ $errors -gt 0 ]]; then
        usage
        exit 1
    fi
    
    # Get absolute paths
    METHFLOW_DIR=$(cd "$METHFLOW_DIR" && pwd)
    mkdir -p "$OUTPUT_DIR"
    OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
    REF=$(realpath "$REF")
}

# Get conda environment path
get_conda_env() {
    local tool=$1
    local conda_yaml="${TOOL_CONFIG[$tool]#*|}"
    conda_yaml="${conda_yaml%%|*}"
    
    local env_check_script="${SCRIPT_DIR}/get_conda_env_path.sh"
    if [[ ! -f "$env_check_script" ]]; then
        log_error "Cannot find environment check script: $env_check_script"
        return 1
    fi
    
    local full_conda_yaml="${METHFLOW_DIR}/${conda_yaml}"
    if [[ ! -f "$full_conda_yaml" ]]; then
        log_error "conda.yaml file does not exist: $full_conda_yaml"
        return 1
    fi
    
    local env_path
    env_path=$("$env_check_script" "$conda_yaml" "$METHFLOW_DIR" 2>/dev/null)
    local exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        local error_msg=$("$env_check_script" "$conda_yaml" "$METHFLOW_DIR" 2>&1 >/dev/null)
        if echo "$error_msg" | grep -q "WARNING"; then
            log_warning "Conda environment has not been created yet"
            log_info "Hint: Run the following command in the container to create the environment:"
            log_info "  cd $METHFLOW_DIR"
            log_info "  snakemake --use-conda --conda-create-envs-only --cores 1"
            if [[ "$SKIP_CONDA_CHECK" == "false" ]]; then
                return 1
            fi
        else
            log_error "Unable to get conda environment path"
            return 1
        fi
    fi

    env_path=$(echo "$env_path" | tr -d '\n\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$env_path" ]]; then
        log_error "Retrieved conda environment path is empty"
        return 1
    fi

    # Check if the environment directory actually exists
    if [[ ! -d "$env_path" ]]; then
        log_error "Conda environment directory does not exist: $env_path"
        log_info "Hint: Run the following command in the container to create the environment:"
        log_info "  cd $METHFLOW_DIR"
        log_info "  snakemake --use-conda --conda-create-envs-only --cores 1"
        if [[ "$SKIP_CONDA_CHECK" == "false" ]]; then
            return 1
        fi
    fi

    log_debug "Conda environment path: $env_path"
    printf '%s' "$env_path"
}

# Set up reference genome (copy or symlink)
setup_reference() {
    local ref=$1
    local target_dir=$2
    local ref_name=$(basename "$ref")
    local target_ref="${target_dir}/${ref_name}"
    
    log_info "Setting up reference genome: $ref_name"
    log_debug "Source file: $ref"
    log_debug "Target: $target_ref"
    
    # Check and handle existing files
    if [[ -e "$target_ref" ]]; then
        if [[ "$COPY_MODE" == "symlink" ]] && [[ -L "$target_ref" ]]; then
            local current_link=$(readlink "$target_ref")
            if [[ "$current_link" == "$ref" ]]; then
                log_debug "Symlink already exists and points to the correct target"
                printf '%s' "$target_ref"
                return 0
            else
                log_warning "Symlink points to a different source file, removing old link..."
                rm -f "$target_ref" "${target_ref}.fai"
            fi
        elif [[ "$COPY_MODE" == "copy" ]] && [[ -f "$target_ref" ]]; then
            log_warning "Target file already exists, reusing existing file"
            printf '%s' "$target_ref"
            return 0
        else
            log_warning "File already exists at target location, removing old file..."
            rm -f "$target_ref" "${target_ref}.fai"
        fi
    fi
    
    # Create copy or symlink
    if [[ "$COPY_MODE" == "copy" ]]; then
        log_info "Copying reference genome to target directory..."
        cp "$ref" "$target_ref"
        [[ -f "${ref}.fai" ]] && cp "${ref}.fai" "${target_ref}.fai"
    else
        log_info "Creating symlink for reference genome..."
        ln -s "$ref" "$target_ref"
        [[ -f "${ref}.fai" ]] && ln -s "${ref}.fai" "${target_ref}.fai"
    fi
    
    log_success "Reference genome has been set up"
    printf '%s' "$target_ref"
}

# Check if index already exists
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
            # Compound path (e.g., Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 or bs/.gem)
            if [[ "$tool" == "gem3" ]]; then
                # gem3 files in bs subdirectory: bs/.gem -> bs/GRCh38.gem
                local dir_part="${suffix%/*}"   # bs
                local file_part="${suffix##*/}" # .gem
                check_file="${output_dir}/${dir_part}/${index_prefix}${file_part}"
            else
                check_file="${output_dir}/${suffix}"
            fi
        elif [[ "$tool" == "whisper" || "$tool" == "gatk" || "$tool" == "bowtie2" || "$tool" == "hisat2" || "$tool" == "hisat-3n" ]]; then
            # These tools use ${index_prefix}${suffix} format (e.g., GRCh38.dict, GRCh38.1.bt2, GRCh38_c2t.3n.CT.1.ht2)
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "gem3" ]]; then
            # gem3 regular index files use ${index_prefix}${suffix} format (e.g., GRCh38.gem)
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "bsbolt" ]]; then
            # bsbolt index files use suffix directly as filename (e.g., BSB_ref.fa, genome_index.pkl)
            check_file="${output_dir}/${suffix}"
        elif [[ "$tool" == "bsmapz" ]]; then
            # bsmapz .fa suffix refers to the reference file itself, not an appended suffix
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
# Generic Template Functions
# =============================================================================

# Generic index building function
build_index_generic() {
    local tool=$1
    local ref=$2
    local output_dir=$3
    local threads=$4
    local conda_env=$5
    
    log_info "Starting to build $tool index..."
    log_info "Reference genome: $(basename "$ref")"
    log_info "Threads: $threads"
    log_debug "Working directory: $output_dir"
    
    cd "$output_dir"
    
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    local half_threads=$((threads / 2))
    [[ $half_threads -lt 1 ]] && half_threads=1
    
    # Get tool configuration
    # Format: env_type|conda_yaml|binary|build_type|args
    local config="${TOOL_CONFIG[$tool]}"
    local env_type="${config%%|*}"
    config="${config#*|}"
    local conda_yaml="${config%%|*}"
    config="${config#*|}"
    local binary="${config%%|*}"
    config="${config#*|}"
    local build_type="${config%%|*}"
    local args="${config#*|}"
    
    # Special handling: hisat-3n requires two build passes
    if [[ "$tool" == "hisat-3n" ]]; then
        build_hisat3n_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # Special handling: gem3 requires two build passes (normal + bisulfite)
    if [[ "$tool" == "gem3" ]]; then
        build_gem3_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # Special handling: batmeth2 requires two build passes (default + RRBS)
    if [[ "$tool" == "batmeth2" ]]; then
        build_batmeth2_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # Special handling: bsmapz does not need index building
    if [[ "$tool" == "bsmapz" ]]; then
        build_bsmapz_ref "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # Special handling: whisper
    if [[ "$tool" == "whisper" ]]; then
        build_whisper_index "$ref" "$output_dir" "$threads"
        return $?
    fi
    
    # Substitute variables in args
    args="${args//\{threads\}/$threads}"
    args="${args//\{half_threads\}/$half_threads}"
    args="${args//\{ref\}/$ref}"
    args="${args//\{index_prefix\}/$index_prefix}"
    args="${args//\{read_length\}/$STROBEALIGN_READ_LENGTH}"
    args="${args//\{dict_file\}/${index_prefix}.dict}"
    
    log_info "Running: $binary $args"

    # Execute command
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
        log_success "$tool index build completed"
        return 0
    else
        log_error "$tool index build failed (exit code: $exit_code)"
        return 1
    fi
}

# Generic index verification function
verify_index_generic() {
    local tool=$1
    local ref=$2
    local output_dir=$3
    
    log_info "Verifying $tool index files..."
    
    local suffixes="${TOOL_SUFFIXES[$tool]}"
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    
    local count=0
    local total=0
    
    for suffix in $suffixes; do
        ((total++))
        local check_file
        if [[ "$suffix" == */* ]]; then
            # Compound path (e.g., Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 or bs/.gem)
            if [[ "$tool" == "gem3" ]]; then
                # gem3 files in bs subdirectory: bs/.gem -> bs/GRCh38.gem
                local dir_part="${suffix%/*}"   # bs
                local file_part="${suffix##*/}" # .gem
                check_file="${output_dir}/${dir_part}/${index_prefix}${file_part}"
            else
                check_file="${output_dir}/${suffix}"
            fi
        elif [[ "$tool" == "whisper" || "$tool" == "gatk" || "$tool" == "bowtie2" || "$tool" == "hisat2" || "$tool" == "hisat-3n" ]]; then
            # These tools use ${index_prefix}${suffix} format (e.g., GRCh38.dict, GRCh38.1.bt2, GRCh38_c2t.3n.CT.1.ht2)
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "gem3" ]]; then
            # gem3 regular index files use ${index_prefix}${suffix} format (e.g., GRCh38.gem)
            check_file="${output_dir}/${index_prefix}${suffix}"
        elif [[ "$tool" == "bsbolt" ]]; then
            # bsbolt index files use suffix directly as filename (e.g., BSB_ref.fa, genome_index.pkl)
            check_file="${output_dir}/${suffix}"
        elif [[ "$tool" == "bsmapz" ]]; then
            # bsmapz .fa suffix refers to the reference file itself, not an appended suffix
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
            log_warning "  $(basename "$check_file") missing"
        fi
    done

    if [[ $count -eq $total ]]; then
        log_success "Index verification passed ($count/$total files exist)"
        return 0
    elif [[ $count -ge $((total - 2)) ]]; then
        log_warning "Index partially exists ($count/$total files exist), may be usable"
        return 0
    else
        log_error "Index verification failed ($count/$total files exist)"
        return 1
    fi
}

# =============================================================================
# Special Tool Handling Functions (requiring multi-step builds)
# =============================================================================

build_hisat3n_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local hisat3n_bin="/opt/MethFlow/resources/hisat-3n/hisat-3n"

    log_info "Starting to build hisat-3n index..."
    log_info "Reference genome: $(basename "$ref")"
    log_info "Threads: $threads"
    log_debug "Working directory: $output_dir"
    log_debug "hisat-3n path: $hisat3n_bin"

    # Check if hisat-3n binary exists
    if [[ ! -f "$hisat3n_bin" ]]; then
        log_error "hisat-3n binary does not exist: $hisat3n_bin"
        return 1
    fi

    # Switch to working directory
    cd "$output_dir"

    local ref_name
    ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"

    # Check if index already exists
    if [[ "$FORCE" == "false" ]] && check_index_exists "hisat-3n" "$ref"; then
        log_warning "hisat-3n index already exists, skipping build (use -f/--force to force rebuild)"
        return 0
    fi

    # Step 1: Build C2T index (--base-change C,T)
    log_info "Step 1/2: Building C2T index..."
    log_info "Running: hisat-3n-build --base-change C,T -p ${threads} $(basename "$ref") ${index_prefix}_c2t"

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
        log_error "hisat-3n C2T index build failed (exit code: $exit_code)"
        return 1
    fi
    log_success "C2T index build completed"

    # Step 2: Build T2C repeat index (--base-change T,C --repeat-index)
    log_info "Step 2/2: Building T2C repeat index..."
    log_info "Running: hisat-3n-build --base-change T,C --repeat-index -p ${threads} $(basename "$ref") ${index_prefix}_t2c"

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
        log_error "hisat-3n T2C repeat index build failed (exit code: $exit_code)"
        return 1
    fi
    log_success "T2C repeat index build completed"

    log_success "All hisat-3n indices built successfully"
    return 0
}

build_gem3_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local gem3_bin="/opt/MethFlow/resources/gem3-mapper/bin/gem-indexer"
    local ref_name=$(basename "$ref")
    local index_prefix="${ref_name%.*}"
    
    log_info "Step 1/2: Building normal index..."
    "$gem3_bin" -i "$ref" -o "${output_dir}/${index_prefix}" -t "$threads" -v 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [gem3]${COLOR_RESET} $line"; done
    
    log_info "Step 2/2: Building bisulfite index..."
    local bs_dir="${output_dir}/bs"
    mkdir -p "$bs_dir"
    ln -sf "$ref" "${bs_dir}/${ref_name}" 2>/dev/null || true
    cd "$bs_dir"
    "$gem3_bin" -i "$ref" -o "${index_prefix}" -b -t "$threads" -v 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [gem3-bs]${COLOR_RESET} $line"; done
    cd "$output_dir"
    
    log_success "gem3 index build completed"
}

build_batmeth2_index() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    local batmeth2_bin="/opt/MethFlow/resources/BatMeth2/bin/BatMeth2"
    
    log_info "Step 1/2: Building default index..."
    "$batmeth2_bin" index -g "$ref" 2>&1 | grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [BatMeth2]${COLOR_RESET} $line"; done
    
    log_info "Step 2/2: Building RRBS index..."
    local rrbs_dir="${output_dir}/rrbs"
    mkdir -p "$rrbs_dir"
    ln -sf "$ref" "${rrbs_dir}/$(basename "$ref")" 2>/dev/null || true
    cd "$rrbs_dir"
    "$batmeth2_bin" index_rrbs -g "$(basename "$ref")" 2>&1 | grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [BatMeth2-RRBS]${COLOR_RESET} $line"; done
    cd "$output_dir"
    
    log_success "BatMeth2 index build completed"
}

build_bsmapz_ref() {
    local ref=$1
    local output_dir=$2
    local threads=$3
    
    log_info "Preparing bsmapz reference genome..."
    local ref_name=$(basename "$ref")
    local target_ref="${output_dir}/${ref_name}"
    
    # Check line length
    local line_len=$(grep -v '^>' "$ref" | grep -v '^$' | head -1 | wc -c)
    if [[ $line_len -le 70 ]]; then
        log_info "Line length meets requirements, creating symlink..."
        ln -sf "$ref" "$target_ref"
    else
        log_info "Line length exceeds 70, reformatting..."
        mamba run -p "$GENOMIC_TOOLS_ENV" seqtk seq -l 70 "$ref" > "$target_ref"
    fi
    
    # Generate .fai
    if [[ ! -f "${target_ref}.fai" ]]; then
        mamba run -p "$GENOMIC_TOOLS_ENV" samtools faidx "$target_ref"
    fi
    
    log_success "bsmapz reference genome preparation completed"
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
    log_info "Running: whisper-index $index_prefix $ref_name $output_dir $temp_dir"
    
    "$whisper_bin" "$index_prefix" "$ref" "$output_dir" "$temp_dir" 2>&1 | \
        grep -iE "(index|build|done|finish)" | while read line; do echo -e "${COLOR_CYAN}    [whisper]${COLOR_RESET} $line"; done
    
    rm -rf "$temp_dir"
    log_success "whisper index build completed"
}

# =============================================================================
# Main Workflow Functions
# =============================================================================

# Process a single tool
process_tool() {
    local tool=$1
    local current=$2
    local total=$3
    
    log_section "Processing tool: $tool"
    log_progress "$current" "$total" "Starting to process $tool"
    
    # Check if tool is supported
    if [[ -z "${TOOL_CONFIG[$tool]:-}" ]]; then
        log_error "Unsupported tool: $tool"
        return 1
    fi
    
    # Parse tool configuration
    local config="${TOOL_CONFIG[$tool]}"
    local env_type="${config%%|*}"
    config="${config#*|}"
    
    # Get conda environment (if needed)
    local conda_env=""
    if [[ "$env_type" == "conda" ]]; then
        conda_env=$(get_conda_env "$tool")
        if [[ $? -ne 0 ]]; then
            log_error "Unable to get conda environment for $tool"
            return 1
        fi
        log_success "Found conda environment"
    else
        log_info "Tool $tool uses container-installed version, skipping conda environment"
    fi
    
    # Create tool output directory
    local tool_dir="${OUTPUT_DIR}/${tool}"
    mkdir -p "$tool_dir"
    log_info "Tool output directory: $tool_dir"
    
    # Set up reference genome
    local target_ref
    target_ref=$(setup_reference "$REF" "$tool_dir")
    
    # Check if index already exists
    if [[ "$FORCE" == "false" ]] && check_index_exists "$tool" "$target_ref"; then
        log_warning "$tool index already exists, skipping build (use -f/--force to force rebuild)"
        verify_index_generic "$tool" "$target_ref" "$tool_dir"
        return 0
    fi
    
    # Build index
    local build_success=false
    build_index_generic "$tool" "$target_ref" "$tool_dir" "$THREADS" "$conda_env"
    
    if [[ $? -eq 0 ]]; then
        verify_index_generic "$tool" "$target_ref" "$tool_dir"
        if [[ $? -eq 0 ]]; then
            build_success=true
        fi
    fi
    
    if [[ "$build_success" == "true" ]]; then
        log_success "$tool index build succeeded"
        return 0
    else
        log_error "$tool index build failed"
        return 1
    fi
}

# Display system information
show_system_info() {
    log_section "System Information"
    if [[ -f /.dockerenv ]] || grep -q docker /proc/1/cgroup 2>/dev/null; then
        log_info "Runtime environment: Docker container"
    else
        log_info "Runtime environment: Host machine"
    fi
    log_info "CPU cores: $(nproc)"
    log_info "Available memory: $(free -h 2>/dev/null | awk '/^Mem:/ {print $2}' || echo "N/A")"
    log_info "Current user: $(whoami)"
    log_info "Working directory: $(pwd)"
}

# =============================================================================
# Main Function
# =============================================================================
main() {
    log_section "Genome Index Building Tool (Template-based Refactored Version)"
    
    # Parse arguments
    parse_args "$@"
    validate_args
    
    # Display system information
    show_system_info
    
    # Display configuration
    log_section "Configuration"
    log_info "Reference genome: $REF"
    log_info "Output directory: $OUTPUT_DIR"
    log_info "Tool list: $TOOLS"
    log_info "Methflow directory: $METHFLOW_DIR"
    log_info "Threads: $THREADS"
    log_info "Copy mode: $COPY_MODE"
    log_info "Force rebuild: $FORCE"
    log_info "strobealign read length: ${STROBEALIGN_READ_LENGTH}bp"
    log_info "Verbose mode: $VERBOSE"
    
    # Split tool list
    IFS=',' read -ra TOOL_ARRAY <<< "$TOOLS"
    local total=${#TOOL_ARRAY[@]}
    log_info "Total tools to process: $total"
    
    # Process each tool
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
    
    # Display summary
    log_section "Build Complete"
    log_info "Total: $total tools"
    log_success "Succeeded: $success_count"
    
    if [[ $fail_count -gt 0 ]]; then
        log_error "Failed: $fail_count"
        exit 1
    fi
    
    log_success "All indices built successfully!"
    log_info "Output directory: $OUTPUT_DIR"
    exit 0
}

# Execute main function
main "$@"