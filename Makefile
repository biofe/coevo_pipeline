.PHONY: install install-dev test lint format clean help

help:  ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package
	pip install -e .

install-dev:  ## Install the package with dev dependencies
	pip install -e ".[dev]"

test:  ## Run all tests
	pytest tests/ -v

test-cov:  ## Run tests with coverage
	pytest tests/ -v --cov=coevo --cov-report=term-missing

lint:  ## Run linter
	python -m flake8 coevo/ tests/ --max-line-length=100

format:  ## Format code with black
	python -m black coevo/ tests/ --line-length=100

clean:  ## Remove build artifacts and caches
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	rm -rf dist/ build/ .coverage htmlcov/

snakemake-dry:  ## Dry-run the Snakemake workflow
	snakemake --dry-run --cores 1

snakemake-run:  ## Run the full Snakemake workflow
	snakemake --cores 8
