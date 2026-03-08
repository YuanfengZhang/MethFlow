#!/usr/bin/env bash
# =============================================================================
# get_conda_env_path.sh
# 根据 conda.yaml 相对路径和 Snakemake 项目根目录，解析 conda 环境目录的绝对地址
# =============================================================================
# 基于 snakemake/deployment/conda.py 中的 Env._get_hash() 方法
# =============================================================================
# Usage: ./get_conda_env_path.sh <conda.yaml 相对路径> <Snakemake 项目根目录>
# Example:
#   ./get_conda_env_path.sh rules/bwa-meth/conda.yaml /mnt/eqa/zhangyuanfeng/methylation/MethFlow
#   ./get_conda_env_path.sh rules/bowtie2/conda.yaml /home/user/project
# =============================================================================

set -euo pipefail

# ---------------------- 函数 ----------------------
usage() {
    cat << EOF
Usage: $0 <conda.yaml 相对路径> <Snakemake 项目根目录>

参数:
  conda.yaml 相对路径    conda.yaml 文件相对于项目根目录的路径
  Snakemake 项目根目录   包含 .snakemake 目录的项目根目录（绝对或相对路径）

示例:
  $0 rules/bwa-meth/conda.yaml /mnt/eqa/zhangyuanfeng/methylation/MethFlow
  $0 rules/bowtie2/conda.yaml /home/user/project
  $0 rules/test/conda.yaml .

输出:
  conda 环境目录的绝对路径

注意:
  - 哈希计算基于 conda.py 中的 Env._get_hash() 方法
  - 包含 .snakemake/conda 目录的 realpath + conda.yaml 内容
  - 自动检测 post-deploy 脚本和 pin 文件
EOF
    exit 1
}

get_conda_env_path() {
    local conda_yaml_rel="$1"
    local project_root="$2"

    # 验证项目根目录
    if [[ ! -d "$project_root" ]]; then
        echo "ERROR: 项目根目录不存在：$project_root" >&2
        return 1
    fi

    # 转换为绝对路径
    local project_root_abs
    project_root_abs="$(cd "$project_root" && pwd)"

    # 验证 conda.yaml 文件
    local conda_yaml_abs="${project_root_abs}/${conda_yaml_rel}"
    if [[ ! -f "$conda_yaml_abs" ]]; then
        echo "ERROR: conda.yaml 文件不存在：$conda_yaml_abs" >&2
        return 1
    fi

    # 验证 .snakemake/conda 目录
    local envs_dir="${project_root_abs}/.snakemake/conda"
    if [[ ! -d "$envs_dir" ]]; then
        echo "ERROR: .snakemake/conda 目录不存在 (请先运行 snakemake 生成环境)" >&2
        return 1
    fi

    # 获取 .snakemake/conda 目录的 realpath（关键！）
    local envs_dir_real
    envs_dir_real="$(readlink -f "$envs_dir")"

    # ---------------------- 计算 MD5 哈希 ----------------------
    # 根据 conda.py _get_hash() 方法：
    # md5hash.update(env_dir.encode())  # 环境目录 realpath
    # md5hash.update(self.content)      # conda.yaml 内容
    # md5hash.update(content_deploy)    # post-deploy 脚本（如果有）
    # md5hash.update(content_pin)       # pin 文件（如果有）

    local hash
    hash=$( {
        # ① 环境目录 realpath
        printf "%s" "$envs_dir_real"
        # ② conda.yaml 原始内容
        cat "$conda_yaml_abs"
        # ③ post-deploy 脚本（如果有）
        local deploy_file="${conda_yaml_abs%.yaml}.post-deploy.sh"
        deploy_file="${deploy_file%.yml}.post-deploy.sh"
        if [[ -f "$deploy_file" ]]; then
            cat "$deploy_file"
        fi
        # ④ pin 文件（如果有）
        local platform
        platform=$(uname -s | tr '[:upper:]' '[:lower:]')
        local pin_file="${conda_yaml_abs%.yaml}.${platform}.pin.txt"
        pin_file="${pin_file%.yml}.${platform}.pin.txt"
        if [[ -f "$pin_file" ]]; then
            cat "$pin_file"
        else
            # 尝试不带平台后缀的 pin 文件
            pin_file="${conda_yaml_abs%.yaml}.pin.txt"
            pin_file="${pin_file%.yml}.pin.txt"
            if [[ -f "$pin_file" ]]; then
                cat "$pin_file"
            fi
        fi
    } | md5sum | cut -d' ' -f1 )

    # ---------------------- 查找环境目录 ----------------------
    # 根据 conda.py address property：
    # hash_candidates = [hash[:8], hash, hash + "_"]

    local hash_short="${hash:0:8}"
    local env_path=""

    for candidate in "${hash}_" "${hash}" "${hash_short}"; do
        if [[ -d "${envs_dir_real}/${candidate}" ]]; then
            env_path="${envs_dir_real}/${candidate}"
            break
        fi
    done

    if [[ -z "$env_path" ]]; then
        # 环境尚未创建，返回预期路径
        echo "WARNING: 环境目录尚未创建，返回预期路径" >&2
        env_path="${envs_dir_real}/${hash_short}_"
    fi

    echo "$env_path"
}

# ---------------------- 主程序 ----------------------
main() {
    if [[ $# -lt 2 ]]; then
        usage
    fi

    local conda_yaml="$1"
    local project_root="$2"

    get_conda_env_path "$conda_yaml" "$project_root"
}

main "$@"