#!/usr/bin/env bash
set -euo pipefail

MODEL_VERSION="1.7-1"
MODEL_NAME="vadr-models-flu-${MODEL_VERSION}"
ARCHIVE="${MODEL_NAME}.tar.gz"
URL="https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu/${MODEL_VERSION}/${ARCHIVE}"
EXPECTED_SHA256="5f09b8d95413251499a2e49a0b93ea119bc96814b4742d92ba55fd3bdadac7ec"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INSTALL_PARENT="${ROOT_DIR}/resources/vadr-models"
INSTALL_DIR="${INSTALL_PARENT}/${MODEL_NAME}"
TMP_DIR="${INSTALL_PARENT}/.tmp-${MODEL_NAME}-$$"
TMP_ARCHIVE="${TMP_DIR}/${ARCHIVE}"

mkdir -p "${INSTALL_PARENT}"

if [[ -f "${INSTALL_DIR}/flu.minfo" && \
      -f "${INSTALL_DIR}/flu.cm" && \
      -f "${INSTALL_DIR}/flu.fa" ]]; then
    echo "VADR influenza models already installed:"
    echo "  ${INSTALL_DIR}"
    exit 0
fi

rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

echo "Downloading VADR influenza models ${MODEL_VERSION}..."
curl -fL "${URL}" -o "${TMP_ARCHIVE}"

if command -v shasum >/dev/null 2>&1; then
    OBSERVED_SHA256="$(shasum -a 256 "${TMP_ARCHIVE}" | awk '{print $1}')"
elif command -v sha256sum >/dev/null 2>&1; then
    OBSERVED_SHA256="$(sha256sum "${TMP_ARCHIVE}" | awk '{print $1}')"
else
    echo "ERROR: neither shasum nor sha256sum is available." >&2
    exit 1
fi

echo "Expected SHA-256: ${EXPECTED_SHA256}"
echo "Observed SHA-256: ${OBSERVED_SHA256}"

if [[ "${OBSERVED_SHA256}" != "${EXPECTED_SHA256}" ]]; then
    echo "ERROR: VADR model archive checksum mismatch." >&2
    exit 1
fi

echo "Checksum verified."

tar -xzf "${TMP_ARCHIVE}" -C "${TMP_DIR}"

EXTRACTED_DIR="${TMP_DIR}/${MODEL_NAME}"

for required in flu.minfo flu.cm flu.fa; do
    if [[ ! -s "${EXTRACTED_DIR}/${required}" ]]; then
        echo "ERROR: required VADR model file missing or empty: ${required}" >&2
        exit 1
    fi
done

rm -rf "${INSTALL_DIR}"
mv "${EXTRACTED_DIR}" "${INSTALL_DIR}"

cat > "${INSTALL_DIR}/WINGS_MODEL_MANIFEST.tsv" <<EOF
component	version	source	archive_sha256
VADR influenza models	${MODEL_VERSION}	${URL}	${EXPECTED_SHA256}
EOF

echo
echo "Installed VADR influenza models:"
echo "  ${INSTALL_DIR}"
echo
echo "VADR model installation complete."