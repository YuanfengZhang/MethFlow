#!/bin/bash

# Setup script for dna_methylation_smk third-party tools
# Usage: ./setup.sh [-c|--choice tool1,tool2,...] [-a|--all] [-h|--help]
# Example: ./setup.sh -c abismal,fame,whisper
#          ./setup.sh --all
#
# NOTE: This script assumes all tools have been cloned as git submodules.
# Run 'git submodule update --init --recursive' before using this script.

set -e

# Base directory for resources
RESOURCES_DIR="$(pwd)/resources"
mkdir -p "$RESOURCES_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Show usage
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Compile third-party tools for dna_methylation_smk
All tools should already be cloned as git submodules.

OPTIONS:
    -c, --choice TOOLS    Comma-separated list of tools to compile
                          Available: abismal,batmeth2,bsgenova,fame,hisat3n,whisper,bioseqzip,rastair
    -a, --all             Compile all tools
    -h, --help            Show this help message

EXAMPLES:
    $0 -c abismal,fame,whisper
    $0 --choice hisat3n
    $0 --all

NOTES:
    - Tools will be compiled in: $RESOURCES_DIR
    - Run 'git submodule update --init --recursive' first to clone all submodules
    - rastair uses pre-compiled binary (downloaded, not compiled)
    - FAME requires manual configuration of CONST.h before compilation
EOF
}

# Install abismal
install_abismal() {
    log_info "Compiling abismal..."
    cd "$RESOURCES_DIR/abismal"
    
    if [ -d "build" ] && [ -f "build/abismal" ]; then
        log_warn "abismal already compiled, skipping..."
        return
    fi
    
    mkdir -p build
    cd build
    ../configure --prefix="$RESOURCES_DIR/abismal"
    make -j$(nproc)
    make install
    log_info "abismal compiled successfully!"
}

# Install BatMeth2
install_batmeth2() {
    log_info "Compiling BatMeth2..."
    cd "$RESOURCES_DIR/BatMeth2"
    
    if [ -f "bin/BatMeth2" ]; then
        log_warn "BatMeth2 already compiled, skipping..."
        return
    fi
    
    # Check dependencies
    if ! command -v gsl-config &> /dev/null; then
        log_error "GSL library not found. Please install libgsl-dev (Ubuntu: apt-get install libgsl-dev)"
        exit 1
    fi
    
    ./configure
    make -j$(nproc)
    log_info "BatMeth2 compiled successfully!"
}

# Install bsgenova
install_bsgenova() {
    log_info "Setting up bsgenova..."
    cd "$RESOURCES_DIR/bsgenova"
    
    # Python scripts, no compilation needed
    # Check if Python dependencies are available
    if ! python3 -c "import numpy" 2>/dev/null; then
        log_warn "numpy not found. bsgenova may require: pip install numpy pandas pysam"
    fi
    
    log_info "bsgenova setup complete! (Python scripts, no compilation needed)"
}

# Install FAME
install_fame() {
    log_info "Compiling FAME..."
    cd "$RESOURCES_DIR/FAME"
    
    if [ -f "FAME" ]; then
        log_warn "FAME already compiled, skipping..."
        return
    fi
    
    log_warn "FAME requires manual configuration of CONST.h before compilation!"
    log_warn "Please edit CONST.h to set CORENUM, CHROMNUM, and READLEN according to your system."
    log_warn "After editing, run: make clean && make"
    
    # Initial compilation with default settings
    make
    log_info "FAME compiled with default settings!"
    log_warn "Remember to recompile after editing CONST.h if needed!"
}

# Install hisat-3n
install_hisat3n() {
    log_info "Compiling hisat-3n..."
    cd "$RESOURCES_DIR/hisat-3n"
    
    if [ -f "hisat-3n" ]; then
        log_warn "hisat-3n already compiled, skipping..."
        return
    fi
    
    make -j$(nproc)
    log_info "hisat-3n compiled successfully!"
}

# Install Whisper
install_whisper() {
    log_info "Compiling Whisper..."
    cd "$RESOURCES_DIR/Whisper/src"
    
    if [ -f "whisper" ]; then
        log_warn "Whisper already compiled, skipping..."
        return
    fi
    
    make -j$(nproc)
    log_info "Whisper compiled successfully!"
}

# Install BioSeqZip
install_bioseqzip() {
    log_info "Compiling BioSeqZip..."
    cd "$RESOURCES_DIR/BioSeqZip"
    
    if [ -d "build" ] && [ -f "build/apps/bioseqzip_collapse" ]; then
        log_warn "BioSeqZip already compiled, skipping..."
        return
    fi
    
    mkdir -p build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    log_info "BioSeqZip compiled successfully!"
}

# Install rastair (pre-compiled binary)
install_rastair() {
    log_info "Downloading rastair pre-compiled binary..."
    cd "$RESOURCES_DIR"
    
    if [ -d "rastair" ]; then
        log_warn "rastair directory already exists, skipping..."
        return
    fi
    
    mkdir -p rastair
    cd rastair
    wget -q https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v2.0.0/rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz
    tar -zxvf rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz
    rm -f rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz
    log_info "rastair installed successfully! (Pre-compiled binary)"
}

# Parse command line arguments
CHOICE=""
INSTALL_ALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--choice)
            CHOICE="$2"
            shift 2
            ;;
        -a|--all)
            INSTALL_ALL=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check if any action specified
if [ "$INSTALL_ALL" = false ] && [ -z "$CHOICE" ]; then
    log_error "No tools specified for installation!"
    show_help
    exit 1
fi

# Check if submodules are initialized
if [ ! -d "$RESOURCES_DIR/abismal" ] || [ ! -d "$RESOURCES_DIR/BatMeth2" ]; then
    log_error "Submodules not found! Please run: git submodule update --init --recursive"
    exit 1
fi

# Install selected tools
if [ "$INSTALL_ALL" = true ]; then
    log_info "Compiling all tools..."
    install_abismal
    install_batmeth2
    install_bsgenova
    install_fame
    install_hisat3n
    install_whisper
    install_bioseqzip
    install_rastair
    log_info "All tools compiled successfully!"
else
    # Parse comma-separated list
    IFS=',' read -ra TOOLS <<< "$CHOICE"
    for tool in "${TOOLS[@]}"; do
        # Trim whitespace
        tool=$(echo "$tool" | xargs)
        case "$tool" in
            abismal)
                install_abismal
                ;;
            batmeth2)
                install_batmeth2
                ;;
            bsgenova)
                install_bsgenova
                ;;
            fame)
                install_fame
                ;;
            hisat3n|hisat-3n)
                install_hisat3n
                ;;
            whisper)
                install_whisper
                ;;
            bioseqzip|BioSeqZip)
                install_bioseqzip
                ;;
            rastair)
                install_rastair
                ;;
            *)
                log_error "Unknown tool: $tool"
                log_info "Available tools: abismal, batmeth2, bsgenova, fame, hisat3n, whisper, bioseqzip, rastair"
                exit 1
                ;;
        esac
    done
    log_info "Selected tools compiled successfully!"
fi

log_info "Setup complete! Tools are in: $RESOURCES_DIR"
