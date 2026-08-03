"""Configuration module for the Mendelian LLM Agent."""
import os
from dotenv import load_dotenv

# --- Path Configuration ---
APP_ROOT = os.path.dirname(os.path.realpath(__file__))
load_dotenv(os.path.join(APP_ROOT, '.env'))

# --- API Keys ---
# Optional, unlike in the upstream working copy. Claude is reached through
# Vertex AI (model id suffixed '@vertex'), which authenticates with
# Application Default Credentials and needs no key at all. Requiring one here
# would make the documented reproduction path impossible. The direct Anthropic
# backend still validates the key when it is actually used.
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

# --- URLs ---
BASE_URL_NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# --- File Paths ---
DATA_DIR = os.path.join(APP_ROOT, 'data')
INPUT_DATA_DIR = os.path.join(DATA_DIR, 'input')
SCORES_FILE = os.path.join(DATA_DIR, 'scores_for_pr_plots.csv')
LOEUF_SMALL_FILE = os.path.join(DATA_DIR, "julia_last_2_small.tsv.bgz")

# Runs directories
RUNS_BASE_DIR = os.path.join(DATA_DIR, 'runs')
MENDELIAN_RUNS_DIR = os.path.join(RUNS_BASE_DIR, 'mendelian')
NOVEL_RUNS_DIR = os.path.join(RUNS_BASE_DIR, 'novel')
PHENOTYPING_RUNS_DIR = os.path.join(RUNS_BASE_DIR, 'phenotyping')

# Separate directories for website and agent runs
WEBSITE_RUNS_DIR = os.path.join(APP_ROOT, 'website_runs')
AGENT_RUNS_DIR = os.path.join(APP_ROOT, 'agent_runs')

# Discovery Potential scores (all_genes_scores.tsv from latest benchmark run)
ALL_GENES_SCORES_FILE = os.path.join(AGENT_RUNS_DIR, 'run_016', 'xgboost', 'fold_5', 'figures', 'all_genes_scores.tsv')

# Legacy (kept for backwards compatibility)
RUNS_DIR = MENDELIAN_RUNS_DIR

# --- API Settings ---
LLM_MAX_RETRIES = 3
LLM_INITIAL_DELAY = 2
API_REQUEST_TIMEOUT = 30

