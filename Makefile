.PHONY: carf-v1 carf-v1-inputs carf-v1-validate carf-v1-validate-all carf-v1-leaderboard carf-v1-sources carf-v1-adapter-smoke carf-v1-download-datasets test manifest verify

carf-v1: carf-v1-inputs carf-v1-validate carf-v1-leaderboard

carf-v1-inputs:
	Rscript --vanilla code/20_prepare_carf_benchmark_v1.R

carf-v1-validate:
	python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/method_scores.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/method_scores.csv
	python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/ground_truth.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/ground_truth.csv
	python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/covariates.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/covariates.csv
	python3 carf_benchmark/scripts/validate_runs.py

carf-v1-validate-all:
	python3 carf_benchmark/scripts/validate_runs.py

carf-v1-leaderboard:
	Rscript --vanilla code/21_run_carf_benchmark_v1.R

carf-v1-sources:
	python3 carf_benchmark/datasets/prepare_public_perturbseq.py materialize-sources --dataset-id all

carf-v1-adapter-smoke:
	python3 carf_benchmark/tests/test_adapters.py

carf-v1-download-datasets:
	python3 carf_benchmark/datasets/prepare_public_perturbseq.py download --dataset-id $${CARF_DATASETS:-all}

test:
	Rscript --vanilla tests/testthat.R
	python3 carf_benchmark/tests/test_adapters.py

manifest:
	Rscript --vanilla code/18_generate_manifest.R

verify:
	Rscript --vanilla code/18_generate_manifest.R --verify
