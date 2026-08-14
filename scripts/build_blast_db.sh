#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

RESOURCE_DIR="${ROOT_DIR}/resources"
ZIP="${RESOURCE_DIR}/fluA_reference.fasta.zip"
DB_DIR="${RESOURCE_DIR}/flu_db"
FASTA="${DB_DIR}/fluA_reference.fasta"
DB_PREFIX="${DB_DIR}/fluA_db"
TMP_DIR="${DB_DIR}/_extract_tmp"
MANIFEST="${DB_DIR}/database_manifest.tsv"

# The reference archive currently comes from the APGAP influenza pipeline release.
GITHUB_REPO="${GITHUB_REPO:-ZooPhy/apgap-influenza-pipeline}"
RELEASE_TAG="${RELEASE_TAG:-v0.1.0}"
ASSET_NAME="${ASSET_NAME:-fluA_reference.fasta.zip}"
ZIP_URL="${ZIP_URL:-https://github.com/${GITHUB_REPO}/releases/download/${RELEASE_TAG}/${ASSET_NAME}}"
BLAST_IMAGE="${BLAST_IMAGE:-ncbi/blast-static:2.17.0}"
BLAST_DB_VERSION="${BLAST_DB_VERSION:-4}"

mkdir -p "${RESOURCE_DIR}" "${DB_DIR}"

sha256_file() {
    local file="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${file}" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "${file}" | awk '{print $1}'
    else
        echo "ERROR: sha256sum or shasum is required to record BLAST database provenance." >&2
        exit 1
    fi
}

if [[ ! -f "${ZIP}" ]]; then
    echo "Compressed reference not found locally: ${ZIP}"
    echo "Attempting to download ${ASSET_NAME} from ${GITHUB_REPO} (${RELEASE_TAG})"

    download_ok=false

    if command -v gh >/dev/null 2>&1; then
        echo "Trying GitHub CLI..."
        if gh release download "${RELEASE_TAG}" \
            --repo "${GITHUB_REPO}" \
            --pattern "${ASSET_NAME}" \
            --dir "${RESOURCE_DIR}" \
            --clobber; then
            [[ -s "${ZIP}" ]] && download_ok=true
        else
            rm -f "${ZIP}"
        fi
    fi

    if [[ "${download_ok}" != true ]]; then
        echo "Trying direct download: ${ZIP_URL}"
        if curl -fL --retry 3 --retry-delay 2 -o "${ZIP}" "${ZIP_URL}"; then
            download_ok=true
        else
            rm -f "${ZIP}"
        fi
    fi

    if [[ "${download_ok}" != true ]]; then
        cat >&2 <<ERROR
ERROR: Unable to obtain ${ASSET_NAME}.

You can fix this by doing one of the following:
  - Authenticate GitHub CLI with: gh auth login
  - Place the archive manually at:
      ${ZIP}
  - Override the source URL:
      ZIP_URL='https://example.org/${ASSET_NAME}' ./scripts/build_blast_db.sh
  - Override the release source:
      GITHUB_REPO='owner/repository' RELEASE_TAG='vX.Y.Z' ./scripts/build_blast_db.sh
ERROR
        exit 1
    fi
fi

if [[ ! -s "${ZIP}" ]]; then
    echo "ERROR: Reference archive is missing or empty: ${ZIP}" >&2
    exit 1
fi

ARCHIVE_SHA256="$(sha256_file "${ZIP}")"

echo "Extracting influenza reference..."
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"
unzip -oq "${ZIP}" -d "${TMP_DIR}"

FASTA_SRC="$(find "${TMP_DIR}" -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \) | head -n 1 || true)"

if [[ -z "${FASTA_SRC}" ]]; then
    echo "ERROR: No FASTA file found inside ${ZIP}" >&2
    rm -rf "${TMP_DIR}"
    exit 1
fi

cp -f "${FASTA_SRC}" "${FASTA}"
FASTA_SHA256="$(sha256_file "${FASTA}")"

# Remove an older database before rebuilding so incompatible database files
# cannot be mixed.
find "${DB_DIR}" -maxdepth 1 -type f \
    \( -name 'fluA_db.n*' -o -name 'fluA_db.*db' \) \
    -delete

BUILD_METHOD=""
MAKEBLASTDB_VERSION=""

build_with_local_blast() {
    echo "Building BLAST database version ${BLAST_DB_VERSION} with local makeblastdb..."
    BUILD_METHOD="local"
    MAKEBLASTDB_VERSION="$(makeblastdb -version 2>&1 | awk '/^makeblastdb:/{print; found=1} END{if(!found) print "unknown"}' | tr '\t' ' ')"
    makeblastdb \
        -in "${FASTA}" \
        -dbtype nucl \
        -out "${DB_PREFIX}"
}

build_with_docker() {
    echo "Building BLAST database version ${BLAST_DB_VERSION} with Docker..."
    BUILD_METHOD="docker"
    MAKEBLASTDB_VERSION="$(docker run --rm "${BLAST_IMAGE}" makeblastdb -version 2>&1 | awk '/^makeblastdb:/{print; found=1} END{if(!found) print "unknown"}' | tr '\t' ' ')"
    docker run --rm \
        -v "${ROOT_DIR}:/workspace" \
        -w /workspace \
        "${BLAST_IMAGE}" \
        makeblastdb \
            -in "/workspace/resources/flu_db/fluA_reference.fasta" \
            -dbtype nucl \
            -out "/workspace/resources/flu_db/fluA_db"
}

build_with_apptainer() {
    echo "Building BLAST database version ${BLAST_DB_VERSION} with Apptainer..."
    BUILD_METHOD="apptainer"
    MAKEBLASTDB_VERSION="$(apptainer exec "docker://${BLAST_IMAGE}" makeblastdb -version 2>&1 | awk '/^makeblastdb:/{print; found=1} END{if(!found) print "unknown"}' | tr '\t' ' ')"
    apptainer exec \
        --bind "${ROOT_DIR}:/workspace" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/workspace/resources/flu_db/fluA_reference.fasta" \
            -dbtype nucl \
            -out "/workspace/resources/flu_db/fluA_db"
}

build_with_singularity() {
    echo "Building BLAST database version ${BLAST_DB_VERSION} with Singularity..."
    BUILD_METHOD="singularity"
    MAKEBLASTDB_VERSION="$(singularity exec "docker://${BLAST_IMAGE}" makeblastdb -version 2>&1 | awk '/^makeblastdb:/{print; found=1} END{if(!found) print "unknown"}' | tr '\t' ' ')"
    singularity exec \
        --bind "${ROOT_DIR}:/workspace" \
        "docker://${BLAST_IMAGE}" \
        makeblastdb \
            -in "/workspace/resources/flu_db/fluA_reference.fasta" \
            -dbtype nucl \
            -out "/workspace/resources/flu_db/fluA_db"
}

if command -v makeblastdb >/dev/null 2>&1; then
    build_with_local_blast
elif command -v docker >/dev/null 2>&1; then
    if ! docker info >/dev/null 2>&1; then
        echo "ERROR: Docker is installed but the Docker daemon is not running." >&2
        echo "Start Docker Desktop and rerun this script." >&2
        rm -rf "${TMP_DIR}"
        exit 1
    fi
    build_with_docker
elif command -v apptainer >/dev/null 2>&1; then
    build_with_apptainer
elif command -v singularity >/dev/null 2>&1; then
    build_with_singularity
else
    echo "ERROR: makeblastdb, Docker, Apptainer, and Singularity were not found." >&2
    rm -rf "${TMP_DIR}"
    exit 1
fi

CREATED_AT_UTC="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

{
    printf 'database_name\tsource_archive\tsource_url\tgithub_repo\trelease_tag\tasset_name\tcreated_at_utc\tarchive_sha256\tfasta_sha256\tmakeblastdb_version\tblast_db_version\tbuild_method\tblast_image\tdb_prefix\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        'fluA_db' \
        "resources/$(basename "${ZIP}")" \
        "${ZIP_URL}" \
        "${GITHUB_REPO}" \
        "${RELEASE_TAG}" \
        "${ASSET_NAME}" \
        "${CREATED_AT_UTC}" \
        "${ARCHIVE_SHA256}" \
        "${FASTA_SHA256}" \
        "${MAKEBLASTDB_VERSION}" \
        "${BLAST_DB_VERSION}" \
        "${BUILD_METHOD}" \
        "${BLAST_IMAGE}" \
        'resources/flu_db/fluA_db'
} > "${MANIFEST}"

rm -rf "${TMP_DIR}"

echo
echo "BLAST database successfully created:"
echo "  ${DB_PREFIX}"
echo "Provenance manifest written:"
echo "  ${MANIFEST}"
echo
echo "Use this value in config.yaml:"
echo "  blast_db: \"resources/flu_db/fluA_db\""
