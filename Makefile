# riposte --- common package tasks. Run `make` (or `make help`) for the list.
# Everything routes through devtools so the whole dependency graph is one tool.
# The slow size/power/exactness simulations run only when NOT_CRAN=true, which the
# relevant targets set for you.

RSCRIPT := Rscript
# Source dir for the sibling fastperm package (Suggests + Remotes). Override to
# the route-b worktree: make fastperm-local FASTPERM_SRC=../fastperm-route-b
FASTPERM_SRC := ../fastperm

.DEFAULT_GOAL := help
.PHONY: help deps document test test-fast check build install fastperm-local vignettes site coverage clean

help: ## List the available targets
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
	  | awk 'BEGIN{FS=":.*?## "}{printf "  make %-14s %s\n", $$1, $$2}'

deps: ## Install development and package dependencies
	$(RSCRIPT) -e 'if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools"); devtools::install_dev_deps(upgrade="never")'

document: ## Regenerate NAMESPACE and the man/*.Rd files from roxygen
	$(RSCRIPT) -e 'devtools::document()'

test: ## Run the full test suite, including the slow simulations
	NOT_CRAN=true $(RSCRIPT) -e 'devtools::test()'

test-fast: ## Run only the fast tests (skip the simulations)
	NOT_CRAN=false $(RSCRIPT) -e 'devtools::test()'

check: document ## R CMD check via devtools (the gate before a change is done)
	$(RSCRIPT) -e 'devtools::check()'

build: document ## Build the source tarball
	$(RSCRIPT) -e 'devtools::build()'

install: document ## Install riposte into the local library
	$(RSCRIPT) -e 'devtools::install(upgrade=FALSE)'

fastperm-local: ## (Re)install sibling fastperm from local source (not GitHub); see FASTPERM_SRC
	$(RSCRIPT) -e 'devtools::install("$(FASTPERM_SRC)", upgrade="never", quick=TRUE)'

vignettes: ## Build the vignettes
	$(RSCRIPT) -e 'devtools::build_vignettes()'

site: ## Build the pkgdown site (needs pkgdown)
	$(RSCRIPT) -e 'pkgdown::build_site()'

coverage: ## Report test coverage (needs covr)
	NOT_CRAN=true $(RSCRIPT) -e 'print(covr::package_coverage())'

clean: ## Remove build artifacts
	rm -rf ..Rcheck *.Rcheck *.tar.gz doc Meta docs
	rm -f src/*.o src/*.so src/*.dll
