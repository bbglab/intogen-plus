#!/bin/bash

# Automatically exit on any failure
set -e

# Those ANSI codes are needed to print with colours
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
GREY='\033[38;5;241m'
RED='\033[0;31m'
# WHITE='\033[0;37m'
BOLD='\033[1m'
RESET='\033[0m' # Reset to default color

# ---------------------------------------------------------------------------------------
# Global variables ----------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# Set CONTAINERS_DIR to the path where this script is located
CONTAINERS_DIR=$(readlink -f "${BASH_SOURCE[0]}" | xargs dirname)

# List of containers (By default all folders with git changes)
CONTAINERS=""

# Whether to run or not docker comands
DRY_RUN="false"

# Whether to skip the build or not
SKIP_BUILD="no"

# Whether to push the containers into the DockerHub or not
PUSH="false"

# What docker binary to use
DOCKER_BIN="docker"

# BuildKit’s advanced features and output capabilities for the Docker builds
export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain

# ---------------------------------------------------------------------------------------
# Function to determine which containers have changes -----------------------------------
# ---------------------------------------------------------------------------------------

changed_containers() {
	( git diff --cached --name-only && git diff --name-only HEAD~1 HEAD ) |
        xargs dirname |
        sed 's|^build/containers/||' | 
        sed -E 's/^(combination|core)/intogen-\1/' |
        sort -u |
        while read -r path; do
            if [ -f "${CONTAINERS_DIR}/$path/Dockerfile" ]; then
                echo "$path"
            fi
        done
}

# ---------------------------------------------------------------------------------------
# Function to determine all the containers available ------------------------------------
# ---------------------------------------------------------------------------------------

all_containers() {
	find -L "${CONTAINERS_DIR}" -name 'Dockerfile' -print0 |
		xargs -0 dirname | xargs basename | sort -u
}

# ---------------------------------------------------------------------------------------
# Function to format a list of containers for display -----------------------------------
# ---------------------------------------------------------------------------------------

format_display() {
	echo "$1" | tr '\n' ',' | sed -E 's/,$//'
}

# ---------------------------------------------------------------------------------------
# Function to format a list of containers into JSON -------------------------------------
# ---------------------------------------------------------------------------------------

format_json() {
	echo "$1" | jq --raw-input . | jq --slurp --compact-output .
}

# ---------------------------------------------------------------------------------------
# Function to build the Docker images ---------------------------------------------------
# ---------------------------------------------------------------------------------------

# This will build each container sequentially.
# It assumes that every container folder will contain a Dockerfile.
# The Docker images will be named after its folder and tagged with the git tag.
# Some folders are just symlinks to upper level folders.
build() {
	container=$1
	tag=$2

	image="${container}:${tag}"
	echo -e "${GREEN}${BOLD}Building Docker image for ${YELLOW}${image}${GREEN} ...${RESET}"
	command="${DOCKER_BIN} build -t ${image} ${CONTAINERS_DIR}/${container}"
	echo -e "${GREY}> ${command}${RESET}"
	if [ "${DRY_RUN}" != "yes" ]; then
        eval "${command}"
    fi
}

# ---------------------------------------------------------------------------------------
# Function to push the containers to the DockerHub --------------------------------------
# ---------------------------------------------------------------------------------------

push() {
	container=$1
	tag=$2

	image="${container}:${tag}"
	echo -e "${GREEN}${BOLD}Pushing Docker image for ${YELLOW}${image}${GREEN} ...${RESET}"
	command="${DOCKER_BIN} push ${image}"
	echo -e "${GREY}> ${command}${RESET}"
	if [ "${DRY_RUN}" != "yes" ]; then
        eval "${command}"
    fi
}

# ---------------------------------------------------------------------------------------
# Parse command line arguments ----------------------------------------------------------
# ---------------------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
	case $1 in
	--changed-containers-as-json)
		format_json "$(changed_containers)"
		exit 0
		;;
	--all-containers-as-json)
		format_json "$(all_containers)"
		exit 0
		;;
	--containers)
        CONTAINERS="${2//,/$'\n'}"
		shift
		shift
		;;
	--all-containers)
		CONTAINERS=$(all_containers)
		shift
		;;
	-n | --dry-run)
		DRY_RUN="yes"
		shift
		;;
	--skip-build)
		SKIP_BUILD="yes"
		shift
		;;
	--push)
		PUSH="yes"
		shift
		;;
	--podman)
		DOCKER_BIN="podman"
		shift
		;;
	*)
		echo -e "${RED}${BOLD}Unknown option: ${YELLOW}$1${RESET}" >&2
		exit 1
		;;
	esac
done

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

if [ -z "${CONTAINERS}" ]; then
	CONTAINERS=$(changed_containers)
fi

CONTAINERS_FMT=$(format_display "$CONTAINERS")

# Set TAG to either the current git tag name or else the git commit short SHA
TAG=$(git describe --tags --exact-match 2>/dev/null || git rev-parse --short HEAD)

echo -e "${CYAN}${BOLD}Containers: ${YELLOW}${CONTAINERS_FMT}${CYAN}${RESET}"

if [ "${SKIP_BUILD}" == "no" ]; then
	for container in ${CONTAINERS}; do
		build "${container}" "${TAG}"
	done

	echo -e "${GREEN}${BOLD}All images built successfully: ${YELLOW}${CONTAINERS_FMT}${CYAN}${RESET}"
fi

if [ "${PUSH}" == "yes" ]; then
	for container in ${CONTAINERS}; do
		push "${container}" "${TAG}"
	done

	echo -e "${GREEN}${BOLD}All images pushed successfully: ${YELLOW}${CONTAINERS_FMT}${CYAN}${RESET}"
fi
