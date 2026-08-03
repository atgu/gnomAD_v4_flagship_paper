"""Data loading utilities for gene scores, GTEx expression and aliases."""
import os
import re
import multiprocessing
import pandas as pd
from config import SCORES_FILE, LOEUF_SMALL_FILE, APP_ROOT, ALL_GENES_SCORES_FILE

GTEX_TPM_FILE = os.path.join(APP_ROOT, "data", "gtex_median_tpm.gct.gz")
ALIASES_FILE = os.path.join(APP_ROOT, "aliases.tsv")
GENCC_FILE = os.path.join(APP_ROOT, "data", "gencc-submissions.tsv")


def _is_quiet() -> bool:
    # Worker processes (ProcessPoolExecutor children) re-import this module
    # and would otherwise spam the parent's tqdm bar with redundant "Loading
    # gene scores..." messages. Auto-silence them. The AGENT_SCORER_QUIET env
    # var can be set explicitly (e.g. via --silent in agent_gene_scorer_v3.py)
    # to also silence the parent process.
    if os.environ.get("AGENT_SCORER_QUIET") == "1":
        return True
    try:
        return multiprocessing.parent_process() is not None
    except Exception:
        return False


def _info(msg: str) -> None:
    if not _is_quiet():
        print(msg)


class DataLoader:
    """Handles loading and caching of gene score data."""
    
    def __init__(self):
        self.gene_scores_df = None
        self.gof_raw_data_df = None
        self.gtex_tpm_df = None
        self.gtex_tissues = None
        self.gencc_df = None
        self.dispo_df = None
        self.alias_map = None  # Maps uppercase alias/previous symbol -> approved symbol
        self.approved_to_aliases = None  # Maps approved symbol uppercase -> set of tokens (including approved)
        self._load_data()
    
    def _load_data(self):
        """Load all data files at initialization."""
        self._load_gene_scores()
        self._load_gof_data()
        self._load_gencc_data()
        self._load_dispo_data()
        # GTEx data is loaded lazily (on first request) due to size
    
    def _load_gene_scores(self):
        """Load gene scores from CSV."""
        _info(f"Loading gene scores from {SCORES_FILE}...")
        try:
            self.gene_scores_df = pd.read_csv(SCORES_FILE)
            _info(f"Loaded {len(self.gene_scores_df)} genes with scores.")
            
            # Check for required columns
            required_cols = ['gene_symbol', 'loeuf_missense_avg', 'obs_missense_avg', 'exp_missense_avg']
            missing_cols = [col for col in required_cols if col not in self.gene_scores_df.columns]
            if missing_cols:
                print(f"WARNING: Missing columns in scores file: {missing_cols}")
        except FileNotFoundError:
            print(f"ERROR: Scores file '{SCORES_FILE}' not found. LOEUF-mis and real_bayes will not be available.")
            self.gene_scores_df = None
        except Exception as e:
            print(f"ERROR: Could not load scores file: {e}")
            self.gene_scores_df = None
    
    def _load_gof_data(self):
        """Pre-load the small GoF data file for performance."""
        _info(f"Pre-loading GoF data from {LOEUF_SMALL_FILE}...")
        try:
            self.gof_raw_data_df = pd.read_csv(LOEUF_SMALL_FILE, sep='\t', compression='gzip')
            _info(f"Successfully pre-loaded {len(self.gof_raw_data_df)} rows for GoF calculations.")
        except FileNotFoundError:
            print(f"WARNING: Small GoF data file not found. GoF calculations will be slower.")
            self.gof_raw_data_df = None
        except Exception as e:
            print(f"ERROR: Could not pre-load small GoF data file: {e}")
            self.gof_raw_data_df = None

    def _load_gencc_data(self):
        """Pre-load GenCC submissions data."""
        _info(f"Loading GenCC data from {GENCC_FILE}...")
        try:
            self.gencc_df = pd.read_csv(GENCC_FILE, sep='\t', dtype=str, low_memory=False)
            # Check for required columns
            required_cols = ['gene_symbol', 'disease_title', 'classification_title']
            missing_cols = [col for col in required_cols if col not in self.gencc_df.columns]
            if missing_cols:
                print(f"WARNING: Missing columns in GenCC file: {missing_cols}")
                self.gencc_df = None
            else:
                # Normalize gene_symbol to uppercase for consistent matching
                self.gencc_df['gene_symbol_upper'] = self.gencc_df['gene_symbol'].str.upper()
                _info(f"Successfully loaded {len(self.gencc_df)} GenCC submissions.")
        except FileNotFoundError:
            print(f"WARNING: GenCC file not found at {GENCC_FILE}. GenCC data will not be available.")
            self.gencc_df = None
        except Exception as e:
            print(f"ERROR: Could not load GenCC data file: {e}")
            self.gencc_df = None

    def _load_dispo_data(self):
        """Load all gene scores from all_genes_scores.tsv and precompute Delta PEPPER."""
        _info(f"Loading gene scores data from {ALL_GENES_SCORES_FILE}...")
        try:
            df = pd.read_csv(ALL_GENES_SCORES_FILE, sep='\t',
                             usecols=['gene_symbol',
                                      'Discovery_Potential', 'Discovery_Potential_pct',
                                      'PEPPER_LLM', 'PEPPER_LLM_pct',
                                      'PEPPER_XGB', 'PEPPER_XGB_pct',
                                      'OMELET_XGB', 'OMELET_XGB_pct'],
                             dtype={'gene_symbol': str})
            df['gene_symbol_upper'] = df['gene_symbol'].str.upper()

            # Precompute Delta PEPPER raw = PEPPER_XGB_pct - PEPPER_LLM_pct
            df['Delta_PEPPER'] = df['PEPPER_XGB_pct'] - df['PEPPER_LLM_pct']
            # Percentile of Delta PEPPER (rank-based, 0-100)
            valid_mask = df['Delta_PEPPER'].notna()
            df.loc[valid_mask, 'Delta_PEPPER_pct'] = (
                df.loc[valid_mask, 'Delta_PEPPER'].rank(pct=True) * 100
            )

            self.dispo_df = df
            _info(f"Loaded gene scores for {len(df)} genes (incl. Delta PEPPER).")
        except FileNotFoundError:
            print(f"WARNING: all_genes_scores.tsv not found at {ALL_GENES_SCORES_FILE}")
            self.dispo_df = None
        except Exception as e:
            print(f"ERROR: Could not load gene scores data: {e}")
            self.dispo_df = None

    def get_discovery_potential(self, gene_name):
        """Return {raw, pct} for a gene, or None if not found."""
        if self.dispo_df is None or not gene_name:
            return None
        try:
            row = self.dispo_df[
                self.dispo_df['gene_symbol_upper'] == gene_name.strip().upper()
            ].iloc[0]
            raw = row['Discovery_Potential']
            pct = row['Discovery_Potential_pct']
            if pd.isna(raw) and pd.isna(pct):
                return None
            return {
                'raw': float(raw) if pd.notna(raw) else None,
                'pct': float(pct) if pd.notna(pct) else None
            }
        except (IndexError, KeyError):
            return None

    def _get_score_pair(self, gene_name, raw_col, pct_col):
        """Generic helper: return {raw, pct} for a given column pair, or None."""
        if self.dispo_df is None or not gene_name:
            return None
        try:
            row = self.dispo_df[
                self.dispo_df['gene_symbol_upper'] == gene_name.strip().upper()
            ].iloc[0]
            raw = row[raw_col]
            pct = row[pct_col]
            if pd.isna(raw) and pd.isna(pct):
                return None
            return {
                'raw': float(raw) if pd.notna(raw) else None,
                'pct': float(pct) if pd.notna(pct) else None
            }
        except (IndexError, KeyError):
            return None

    def get_pepper_llm(self, gene_name):
        return self._get_score_pair(gene_name, 'PEPPER_LLM', 'PEPPER_LLM_pct')

    def get_pepper_xgb(self, gene_name):
        return self._get_score_pair(gene_name, 'PEPPER_XGB', 'PEPPER_XGB_pct')

    def get_omelet_xgb(self, gene_name):
        return self._get_score_pair(gene_name, 'OMELET_XGB', 'OMELET_XGB_pct')

    def get_delta_pepper(self, gene_name):
        return self._get_score_pair(gene_name, 'Delta_PEPPER', 'Delta_PEPPER_pct')

    def _load_alias_map(self):
        """Lazy-load HGNC aliases TSV and build alias->approved map (first hit wins)."""
        if self.alias_map is not None:
            return

        if not os.path.exists(ALIASES_FILE):
            print(f"WARNING: aliases file not found at {ALIASES_FILE}")
            self.alias_map = {}
            return

        try:
            df = pd.read_csv(ALIASES_FILE, sep='\t', dtype=str).fillna('')
        except Exception as e:
            print(f"ERROR: Could not load aliases file {ALIASES_FILE}: {e}")
            self.alias_map = {}
            return

        alias_map = {}
        approved_to_aliases = {}

        def add_alias(token, approved_symbol):
            alias_upper = token.strip().upper()
            if not alias_upper:
                return
            if alias_upper not in alias_map:
                alias_map[alias_upper] = approved_symbol
            approved_upper = approved_symbol.upper()
            approved_to_aliases.setdefault(approved_upper, set()).add(alias_upper)

        for _, row in df.iterrows():
            approved = (row.get('Approved symbol') or '').strip()
            if not approved:
                continue
            approved_symbol = approved
            add_alias(approved_symbol, approved_symbol)

            prev_symbols = (row.get('Previous symbols') or '')
            alias_symbols = (row.get('Alias symbols') or '')

            for cell in (prev_symbols, alias_symbols):
                if not cell:
                    continue
                for token in re.split(r'[,;]', cell):
                    add_alias(token, approved_symbol)

        self.alias_map = alias_map
        self.approved_to_aliases = approved_to_aliases

    def _resolve_gene_symbol(self, gene_name):
        """Return approved symbol if gene_name matches an alias/previous symbol."""
        if not gene_name:
            return gene_name
        self._load_alias_map()
        return self.alias_map.get(gene_name.strip().upper(), gene_name)

    def _deduplicate_gencc_entries(self, entries):
        """
        Deduplicate GenCC entries by name inclusion.
        - If one disease name contains another (or is equal), merge them
        - Keep the shortest name
        - Merge and deduplicate classification labels (keep unique ones)
        
        Args:
            entries: List of dicts with 'disease_title' and 'classification_title'
            
        Returns:
            List of deduplicated entries
        """
        if not entries:
            return []
        
        # Normalize names for comparison (case-insensitive)
        def normalize_name(name):
            return name.lower().strip()
        
        # Group entries by matching names (transitive closure)
        merged_groups = []
        processed_indices = set()
        
        for i, entry1 in enumerate(entries):
            if i in processed_indices:
                continue
            
            # Start a new group with this entry
            group = {
                'names': [entry1['disease_title']],
                'classifications': {entry1['classification_title']}  # Use set for uniqueness
            }
            processed_indices.add(i)
            
            # Find all entries that match this one (transitive)
            changed = True
            while changed:
                changed = False
                for j, entry2 in enumerate(entries):
                    if j in processed_indices:
                        continue
                    
                    name1_normalized = normalize_name(entry2['disease_title'])
                    
                    # Check if entry2's name matches any name in the current group
                    for existing_name in group['names']:
                        existing_normalized = normalize_name(existing_name)
                        
                        # Check if one name contains the other (or is equal)
                        if (name1_normalized == existing_normalized or
                            name1_normalized in existing_normalized or
                            existing_normalized in name1_normalized):
                            # Match found - add to group
                            group['names'].append(entry2['disease_title'])
                            group['classifications'].add(entry2['classification_title'])
                            processed_indices.add(j)
                            changed = True
                            break
            
            merged_groups.append(group)
        
        # Build final result list
        result_list = []
        for group in merged_groups:
            # Keep the shortest name
            shortest_name = min(group['names'], key=len)
            
            # Get unique classifications (sorted for consistency)
            unique_classifications = sorted(list(group['classifications']))
            
            # Create a single entry with all classifications
            result_list.append({
                'disease_title': shortest_name,
                'classification_title': unique_classifications[0] if unique_classifications else 'N/A',
                'all_classifications': unique_classifications  # Store all classifications
            })
        
        return result_list
    
    def get_aliases_for_gene(self, gene_name):
        """
        Return list of aliases/previous symbols for the approved gene (excludes approved symbol).
        """
        if not gene_name:
            return []
        self._load_alias_map()
        resolved = self._resolve_gene_symbol(gene_name)
        approved_upper = (resolved or "").strip().upper()
        if not approved_upper:
            return []
        tokens = list(self.approved_to_aliases.get(approved_upper, set()))
        # Remove the approved symbol itself (keep only alternatives)
        aliases = [t for t in tokens if t != approved_upper]
        return aliases
    
    def get_gene_data(self, gene_name):
        """
        Retrieve all relevant scores (score, obs, exp) for a gene for all LOEUF versions.
        """
        if self.gene_scores_df is None:
            return None

        try:
            gene_info = self.gene_scores_df[self.gene_scores_df['gene_symbol'].str.upper() == gene_name.upper()].iloc[0]
            
            all_scores = {
                "loeuf_mis": {
                    "score": gene_info.get('loeuf_missense_avg'),
                    "obs": gene_info.get('obs_missense_avg'),
                    "exp": gene_info.get('exp_missense_avg')
                },
                "loeuf_v2": {
                    "score": gene_info.get('oe_lof_upper_v2'),
                    "obs": gene_info.get('obs_lof_v2'),
                    "exp": gene_info.get('exp_lof_v2')
                },
                "loeuf_v4": {
                    "score": gene_info.get('loeuf_linear_new_loftee_99_5_adj_r'),
                    "obs": gene_info.get('linear__new_loftee_99_5__adj_r_obs'),
                    "exp": gene_info.get('linear__new_loftee_99_5__adj_r_exp')
                }
            }
            return all_scores
        except (IndexError, KeyError):
            return None
    
    def _load_gtex_data(self):
        """Lazy load GTEx TPM data."""
        if self.gtex_tpm_df is not None:
            return
        
        _info(f"Loading GTEx TPM data from {GTEX_TPM_FILE}...")
        try:
            # Skip first 2 lines (header info), read the data
            self.gtex_tpm_df = pd.read_csv(GTEX_TPM_FILE, sep='\t', skiprows=2, compression='gzip')
            # Store tissue names (all columns except Name and Description)
            self.gtex_tissues = [col for col in self.gtex_tpm_df.columns if col not in ['Name', 'Description']]
            _info(f"Loaded GTEx data for {len(self.gtex_tpm_df)} genes, {len(self.gtex_tissues)} tissues.")
        except FileNotFoundError:
            print(f"WARNING: GTEx TPM file not found at {GTEX_TPM_FILE}")
            self.gtex_tpm_df = pd.DataFrame()
            self.gtex_tissues = []
        except Exception as e:
            print(f"ERROR: Could not load GTEx data: {e}")
            self.gtex_tpm_df = pd.DataFrame()
            self.gtex_tissues = []
    
    def get_top_tissues(self, gene_name, top_n=10):
        """
        Get top N tissues by expression for a gene.
        
        Returns:
            list of dicts: [{tissue: str, tpm: float}, ...] sorted by TPM descending
        """
        self._load_gtex_data()
        
        if self.gtex_tpm_df is None or len(self.gtex_tpm_df) == 0:
            return []
        
        # Build search candidates: approved + any aliases pointing to it + raw input
        search_symbols = []
        seen = set()
        resolved = self._resolve_gene_symbol(gene_name)
        if resolved:
            for candidate in self.approved_to_aliases.get(resolved.upper(), {resolved}):
                if candidate not in seen:
                    search_symbols.append(candidate)
                    seen.add(candidate)
        if gene_name:
            raw_upper = gene_name.strip().upper()
            if raw_upper not in seen:
                search_symbols.append(gene_name)
                seen.add(raw_upper)

        try:
            gene_row = None
            for symbol in search_symbols:
                candidate = self.gtex_tpm_df[
                    self.gtex_tpm_df['Description'].str.upper() == symbol.upper()
                ]
                if not candidate.empty:
                    gene_row = candidate.iloc[0]
                    break

            if gene_row is None:
                return []
            
            # Get all tissue values
            tissue_values = []
            for tissue in self.gtex_tissues:
                tpm = gene_row.get(tissue, 0)
                if pd.notna(tpm):
                    tissue_values.append({'tissue': tissue, 'tpm': float(tpm)})
            
            # Sort by TPM descending and take top N
            tissue_values.sort(key=lambda x: x['tpm'], reverse=True)
            return tissue_values[:top_n]
            
        except Exception as e:
            print(f"ERROR: Could not get top tissues for {gene_name}: {e}")
            return []

    def get_gencc_diseases(self, gene_name):
        """
        Get list of unique disease entries from GenCC for a given gene.
        
        Args:
            gene_name (str): Gene symbol to search for
            
        Returns:
            list: List of dicts with 'disease_title' and 'classification_title', 
                  empty list if gene not found or data unavailable
        """
        if self.gencc_df is None or len(self.gencc_df) == 0:
            return []
        
        if not gene_name:
            return []
        
        try:
            # Search by gene_symbol (case-insensitive)
            gene_upper = gene_name.strip().upper()
            matching_rows = self.gencc_df[
                self.gencc_df['gene_symbol_upper'] == gene_upper
            ]
            
            if matching_rows.empty:
                return []
            
            # Extract disease_title and classification_title, filtering out empty/NaN values
            raw_entries = []
            
            for _, row in matching_rows.iterrows():
                disease_title = str(row.get('disease_title', '')).strip()
                classification_title = str(row.get('classification_title', '')).strip()
                
                # Remove quotes if present (pandas might preserve quotes from TSV)
                disease_title = disease_title.strip('"\'')
                classification_title = classification_title.strip('"\'')
                
                # Skip if disease_title is empty
                if not disease_title or disease_title == 'nan' or disease_title == '':
                    continue
                
                # Normalize classification_title - use 'N/A' if empty or invalid
                if not classification_title or classification_title == 'nan' or classification_title == '':
                    classification_title = 'N/A'
                
                raw_entries.append({
                    'disease_title': disease_title,
                    'classification_title': classification_title
                })
            
            # Deduplicate entries by name inclusion and merge classifications
            result_list = self._deduplicate_gencc_entries(raw_entries)
            
            # Sort by evidence level (classification_title) first, then by disease_title
            # Evidence order from strongest to weakest (using exact values from GenCC file):
            evidence_order = {
                'Definitive': 1,
                'Strong': 2,
                'Moderate': 3,
                'Supportive': 4,
                'Limited': 5,
                'Disputed Evidence': 6,
                'Refuted Evidence': 7,
                'Animal': 8,
                'No Known Disease Relationship': 9
            }
            
            def sort_key(item):
                classification = item['classification_title']
                # Normalize classification for matching (case-insensitive, strip whitespace)
                classification_normalized = classification.strip() if classification else ''
                evidence_rank = evidence_order.get(classification_normalized, 99)  # Unknown classifications go last
                disease_name = item['disease_title']
                return (evidence_rank, disease_name)
            
            result_list.sort(key=sort_key)
            
            return result_list
            
        except Exception as e:
            print(f"ERROR: Could not get GenCC diseases for {gene_name}: {e}")
            return []


# Global data loader instance
data_loader = DataLoader()

