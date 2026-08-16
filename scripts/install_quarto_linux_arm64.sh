#!/usr/bin/env bash
set -euo pipefail

QUARTO_VERSION="1.9.38"
QUARTO_SHA256="75fbc5c1121ffe65e564e9d24711db2ad8f617f9552f5dc7d8a06307d72dde38"

INSTALL_ROOT="${1:-software}"
INSTALL_DIR="${INSTALL_ROOT}/quarto-${QUARTO_VERSION}"
ARCHIVE="quarto-${QUARTO_VERSION}-linux-arm64.tar.gz"
URL="https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/${ARCHIVE}"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "ERROR: This installer is for Linux ARM64 only." >&2
    exit 1
fi

case "$(uname -m)" in
    aarch64|arm64)
        ;;
    *)
        echo "ERROR: Unsupported architecture: $(uname -m). Expected aarch64 or arm64." >&2
        exit 1
        ;;
esac

if [[ -x "${INSTALL_DIR}/bin/quarto" ]]; then
    echo "Quarto ${QUARTO_VERSION} is already installed:"
    "${INSTALL_DIR}/bin/quarto" --version
    exit 0
fi

mkdir -p "${INSTALL_ROOT}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

echo "Downloading Quarto ${QUARTO_VERSION} for Linux ARM64..."
curl -fL "${URL}" -o "${tmpdir}/${ARCHIVE}"

echo "${QUARTO_SHA256}  ${tmpdir}/${ARCHIVE}" | sha256sum -c -

mkdir -p "${INSTALL_DIR}"
tar -xzf "${tmpdir}/${ARCHIVE}" -C "${INSTALL_DIR}" --strip-components=1

if [[ ! -x "${INSTALL_DIR}/bin/quarto" ]]; then
    echo "ERROR: Quarto executable was not installed correctly." >&2
    exit 1
fi

echo
echo "Installed Quarto ${QUARTO_VERSION}:"
"${INSTALL_DIR}/bin/quarto" --version

echo
echo "Add this directory to PATH before running WINGS:"
printf '  export PATH="%s/bin:$PATH"\n' "$(cd "${INSTALL_DIR}" && pwd)"
