CURRENT_DIR := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))

BUILD := $(CURRENT_DIR)/build.sh
BUILD_ARGS ?=

SHELLCHECK_BIN ?= shellcheck
SHELLCHECK_ARGS ?= 

HADOLINT_IGNORES := DL3059 # Multiple consecutive `RUN` instructions
HADOLINT_IGNORES += DL3003 # Use WORKDIR to switch to a directory
HADOLINT_IGNORES += DL3042 # Avoid use of cache directory with pip

# HADOLINT_BIN := docker run --rm -i hadolint/hadolint hadolint
HADOLINT_BIN ?= hadolint
HADOLINT_ARGS ?= $(foreach ignore,$(HADOLINT_IGNORES),--ignore $(ignore))

# Those ANSI codes are needed to print with colours
GREEN := \033[0;32m
WHITE := \033[0;37m
BOLD := \033[1m
BOLD_GREEN := $(GREEN)$(BOLD)
BOLD_WHITE := $(WHITE)$(BOLD)
RESET := \033[0m

.DEFAULT_GOAL = checks

.PHONY: checks
checks: shellcheck hadolint

.PHONY: shellcheck
shellcheck:
	@echo "$(BOLD_GREEN)Running cheks for shell scripts ...$(RESET)"
	set -e; \
	find -L $(CURRENT_DIR) -name '*.sh' | sed 's|$(CURRENT_DIR)/||' \
		| while read -r script; do \
		  echo "$(BOLD_WHITE)=> $${script}$(RESET)"; \
		  $(SHELLCHECK_BIN) $(SHELLCHECK_ARGS) $${script}; \
		  done

.PHONY: hadolint
hadolint:
	@echo "$(BOLD_GREEN)Running checks for Dockerfile ...$(RESET)"
	set -e; \
	find -L $(CURRENT_DIR) -name 'Dockerfile' | sed 's|$(CURRENT_DIR)/||' \
		| while read -r dockerfile; do \
		  echo "$(BOLD_WHITE)=> $${dockerfile}$(RESET)"; \
		  $(HADOLINT_BIN) $(HADOLINT_ARGS) - < $${dockerfile}; \
		  done

.PHONY: build-all
build-all:
	$(BUILD) --all-containers $(BUILD_ARGS)

.PHONY: build-changes
build-changes:
	$(BUILD) $(BUILD_ARGS)

.PHONY: push-all
push-all:
	$(BUILD) --push --all-containers $(BUILD_ARGS)

.PHONY: push-changes
push-changes:
	$(BUILD) --push $(BUILD_ARGS)
