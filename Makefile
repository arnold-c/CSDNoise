include config.mk

.PHONY: all
all: ensemble-targets tycho-targets test-targets

# Ensemble targets
ENSEMBLE_TARGETS = ensemble-sim ensemble-sim-ews-visualization
.PHONY: $(ENSEMBLE_TARGETS) ensemble-targets
$(ENSEMBLE_TARGETS): %: tmp/%
ensemble-targets: $(ENSEMBLE_TARGETS)

tmp/ensemble-sim: scripts/ensemble-sim.jl
	julia $^
	@touch $@

tmp/ensemble-sim-ews-visualization: scripts/ensemble-sim_ews-visualization.jl tmp/ensemble-sim
	julia $<
	@touch $@

# Tycho visualization targets
TYCHO_TARGETS = tycho-visualization tycho-brett-visualization-R tycho-noise-R tycho-data-prep
.PHONY: $(TYCHO_TARGETS) tycho-targets
$(TYCHO_TARGETS): %: tmp/%
tmp/tycho-visualization: scripts/tycho-visualization.jl tmp/tycho-data-prep
	julia $<
	@touch $@

tmp/tycho-data-prep: scripts/tycho-cleaning.R
	rig switch $(R_VERSION)
	Rscript $<
	@touch $@

tmp/tycho-brett-visualization-R: scripts/tycho-brett-visualization.R tmp/tycho-data-prep
	rig switch $(R_VERSION)
	Rscript $^
	@touch $@

tmp/tycho-noise-R: scripts/tycho-noise.R tmp/tycho-brett-visualization-R
	rig switch $(R_VERSION)
	Rscript $^
	@touch $@

# Test targets
TEST_TARGETS = runtests
.PHONY: $(TEST_TARGETS) test-targets
$(TEST_TARGETS): %: tmp/%
test-targets: $(TEST_TARGETS)

tmp/runtests: test/runtests.jl ensemble-targets
	julia $<
	@touch $@

# Cleaning targets
.PHONY: clean-tests clean-all clean-tmp clean-plots clean-all-ensemble clean-ensemble-sims clean-ensemble-sims-ews-visualization clean-tycho-visualization clean-tycho-data-prep clean-tycho-brett-visualization
clean-all: clean-tests clean-tmp clean-plots clean-all-ensemble clean-ensemble-sims-ews-visualization clean-tycho-visualization clean-tycho-data-prep clean-tycho-brett-visualization

clean-tests:
	@echo "cleaning tests"
	$(shell fd -g 'test' 'tmp/' | xargs rm)

clean-tmp:
	@echo "cleaning all tmp files"
	$(shell rm -rf tmp)

clean-plots:
	@echo "cleaning plot output files"
	$(shell fd -g '*.png' 'plots/' -HI | xargs rm -r)

clean-all-ensemble: clean-ensemble-sims

clean-ensemble-sims:
	@echo "cleaning all ensemble output files"
	$(shell fd . 'out/ensemble/' -td -HI | xargs rm -r)
	@echo "cleaning ensemble simulation tmp files"
	$(shell fd -g 'ensemble-sim' 'tmp/' | xargs rm)

clean-ensemble-sim-ews-visualization:
	@echo "cleaning ensemble EWS visualization plots"
	$(shell fd 'burnin-time' 'plots/ensemble/' | xargs rm -r)
	@echo "cleaning ensemble EWS visualization tmp files"
	$(shell fd -g 'ensemble-sim-ews-visualization' 'tmp/' | xargs rm)

clean-tycho-visualization:
	@echo "cleaning tycho visualization"
	$(shell fd -g 'tycho' 'plots/' | xargs rm -r)

clean-tycho-data-prep:
	@echo "cleaning tycho data prep"
	$(shell fd . 'out/tycho/R-scripts/' | xargs rm)
	$(shell fd 'tycho.*.csv' 'out/tycho/' | xargs rm)
	@echo "cleaning tycho data prep tmp files"
	$(shell fd -g 'tycho-data-prep' 'tmp/' | xargs rm)

clean-tycho-brett-visualization:
	@echo "cleaning tycho Brett visualization"
	$(shell fd -g 'tycho-brett' 'plots/' | xargs rm -r)
