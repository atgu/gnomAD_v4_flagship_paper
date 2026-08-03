"""Service for managing analysis run persistence."""
import os
import json
import glob
import re


class RunPersistenceService:
    """Service to handle saving and loading analysis runs."""
    
    def __init__(self, runs_dir):
        """
        Initialize the persistence service.
        
        Args:
            runs_dir: Directory where runs should be saved
        """
        self.runs_dir = runs_dir
        # Ensure directory exists
        if not os.path.exists(self.runs_dir):
            os.makedirs(self.runs_dir)
    
    def save_mendelian_run(self, gene_name, model, num_papers, keywords, 
                          prompt_template, articles, llm_response_text, parsed_response):
        """
        Save a Mendelian analysis run to JSON file.
        
        Args:
            gene_name: Name of the gene analyzed
            model: LLM model used
            num_papers: Number of papers retrieved
            keywords: Search keywords used
            prompt_template: Prompt template used
            articles: PubMed articles retrieved
            llm_response_text: Raw LLM response
            parsed_response: Parsed response with all results
            
        Returns:
            Filepath where the run was saved, or None if save failed
        """
        try:
            # Get next run number
            run_num = self.get_next_run_number(gene_name)
            
            filename = f"{gene_name}_run{run_num:02d}.json"
            filepath = os.path.join(self.runs_dir, filename)

            run_data = {
                "gene_symbol": gene_name,
                "run_parameters": {
                    "model": model,
                    "num_papers": num_papers,
                    "keywords": keywords,
                    "prompt": prompt_template
                },
                "pubmed_articles_sent_to_llm": articles,
                "llm_raw_response": llm_response_text,
                "final_output": parsed_response
            }

            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(run_data, f, indent=4)
            
            print(f"INFO: Saved detailed run output to '{filepath}'")
            return filepath

        except Exception as e:
            print(f"ERROR: Could not save run output JSON for {gene_name}. Error: {e}")
            return None
    
    def save_novel_run(self, gene_name, phenotype, synonyms, model,
                      articles, llm_response, parsed_result,
                      prompt_templates=None, prompt_texts=None,
                      parameters=None, api_response=None,
                      pubmed_pmids=None):
        """
        Save a novel association analysis run to JSON file.
        
        Args:
            gene_name: Name of the gene
            phenotype: Phenotype analyzed
            synonyms: Phenotype synonyms used
            model: LLM model used
            articles: PubMed articles retrieved
            llm_response: Raw LLM response
            parsed_result: Parsed result
            
        Returns:
            Filepath where the run was saved, or None if save failed
        """
        try:
            run_num = 1
            safe_phenotype = phenotype.replace(" ", "_").replace("/", "_")
            
            # Find next available run number
            while os.path.exists(os.path.join(self.runs_dir, f"{gene_name}_{safe_phenotype}_novel_run{run_num:02d}.json")):
                run_num += 1
            
            filename = f"{gene_name}_{safe_phenotype}_novel_run{run_num:02d}.json"
            filepath = os.path.join(self.runs_dir, filename)
            
            run_data = {
                "gene_name": gene_name,
                "phenotype": phenotype,
                "synonyms": synonyms,
                "run_parameters": {
                    "model": model,
                    "temperature": (parameters or {}).get("temperature"),
                    "batch_size": (parameters or {}).get("batch_size")
                },
                "prompt_templates": prompt_templates or {},
                "prompt_texts": prompt_texts or {},
                "pubmed_articles": articles,
                "pubmed_pmids": pubmed_pmids or [],
                "llm_raw_response": llm_response,
                "final_output": parsed_result,
                "api_response": api_response
            }
            
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(run_data, f, indent=4)
            
            print(f"INFO: Saved novel association run output to '{filepath}'")
            return filepath
            
        except Exception as e:
            print(f"ERROR: Could not save novel run output JSON. Error: {e}")
            return None
    
    def get_latest_run_number(self):
        """
        Get the highest run number in the runs directory.
        
        Returns:
            int: Latest run number, or 0 if no runs exist
        """
        pattern = re.compile(r'^run_(\d+)$')
        run_numbers = []
        
        try:
            for item in os.listdir(self.runs_dir):
                item_path = os.path.join(self.runs_dir, item)
                if os.path.isdir(item_path):
                    match = pattern.match(item)
                    if match:
                        run_numbers.append(int(match.group(1)))
        except Exception as e:
            print(f"ERROR: Could not list runs: {e}")
            return 0
        
        return max(run_numbers) if run_numbers else 0
    
    def get_latest_run_path(self):
        """
        Get the path to the latest run directory.
        
        Returns:
            str: Path to latest run, or None if no runs exist
        """
        latest_num = self.get_latest_run_number()
        if latest_num == 0:
            return None
        return os.path.join(self.runs_dir, f"run_{latest_num:03d}")
    
    def create_new_run(self, prompts, config):
        """
        Create a new run directory with prompts and config.
        
        Args:
            prompts: Dictionary of all prompts
            config: Dictionary of configuration parameters
            
        Returns:
            str: Path to the new run directory
        """
        next_num = self.get_latest_run_number() + 1
        run_path = os.path.join(self.runs_dir, f"run_{next_num:03d}")
        
        try:
            # Create run directory and results subdirectory
            os.makedirs(run_path, exist_ok=True)
            os.makedirs(os.path.join(run_path, "results"), exist_ok=True)
            
            # Save prompts
            prompts_path = os.path.join(run_path, "prompts.json")
            with open(prompts_path, 'w', encoding='utf-8') as f:
                json.dump(prompts, f, indent=2)
            
            # Save config
            config_path = os.path.join(run_path, "config.json")
            with open(config_path, 'w', encoding='utf-8') as f:
                json.dump(config, f, indent=2)
            
            print(f"INFO: Created new run: {run_path}")
            return run_path
            
        except Exception as e:
            print(f"ERROR: Could not create new run: {e}")
            return None
    
    def load_run_prompts(self, run_path):
        """
        Load prompts from a run directory.
        
        Args:
            run_path: Path to run directory
            
        Returns:
            dict: Dictionary of prompts, or None if loading failed
        """
        try:
            prompts_path = os.path.join(run_path, "prompts.json")
            if not os.path.exists(prompts_path):
                return None
            
            with open(prompts_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"ERROR: Could not load prompts from {run_path}: {e}")
            return None
    
    def load_run_config(self, run_path):
        """
        Load config from a run directory.
        
        Args:
            run_path: Path to run directory
            
        Returns:
            dict: Dictionary of config, or None if loading failed
        """
        try:
            config_path = os.path.join(run_path, "config.json")
            if not os.path.exists(config_path):
                return None
            
            with open(config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"ERROR: Could not load config from {run_path}: {e}")
            return None
    
    def save_gene_result(self, run_path, gene_name, result_data):
        """
        Save a gene result to a run's results directory.
        
        Args:
            run_path: Path to run directory
            gene_name: Name of the gene
            result_data: Dictionary with analysis results
            
        Returns:
            str: Path to saved file, or None if save failed
        """
        try:
            results_dir = os.path.join(run_path, "results")
            os.makedirs(results_dir, exist_ok=True)
            
            # Sanitize gene name for filename
            safe_gene_name = "".join(c for c in gene_name if c.isalnum() or c in ('-', '_')).strip()
            if not safe_gene_name:
                safe_gene_name = "unknown_gene"
            
            filepath = os.path.join(results_dir, f"{safe_gene_name}.json")
            
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(result_data, f, indent=2)
            
            return filepath
            
        except Exception as e:
            print(f"ERROR: Could not save result for {gene_name}: {e}")
            return None

    def save_gene_raw_llm(self, run_path, gene_name, raw_llm_data):
        """
        Save raw LLM artifacts (prompts + raw responses) for a gene alongside results.

        We intentionally keep this separate from the main results JSON to avoid
        bloating the primary payload consumed by the web UI.

        File: <run_path>/results/<GENE>.raw_llm.json
        """
        try:
            results_dir = os.path.join(run_path, "results")
            os.makedirs(results_dir, exist_ok=True)

            safe_gene_name = "".join(c for c in gene_name if c.isalnum() or c in ('-', '_')).strip()
            if not safe_gene_name:
                safe_gene_name = "unknown_gene"

            filepath = os.path.join(results_dir, f"{safe_gene_name}.raw_llm.json")

            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(raw_llm_data, f, indent=2)

            return filepath
        except Exception as e:
            print(f"ERROR: Could not save raw LLM data for {gene_name}: {e}")
            return None
    
    def load_gene_result(self, run_path, gene_name):
        """
        Load a gene result from a run's results directory.
        
        Args:
            run_path: Path to run directory
            gene_name: Name of the gene
            
        Returns:
            dict: Gene result data, or None if not found
        """
        try:
            # Sanitize gene name for filename
            safe_gene_name = "".join(c for c in gene_name if c.isalnum() or c in ('-', '_')).strip()
            if not safe_gene_name:
                safe_gene_name = "unknown_gene"
            
            filepath = os.path.join(run_path, "results", f"{safe_gene_name}.json")
            
            if not os.path.exists(filepath):
                return None
            
            with open(filepath, 'r', encoding='utf-8') as f:
                return json.load(f)
                
        except Exception as e:
            print(f"ERROR: Could not load result for {gene_name}: {e}")
            return None
    
    def config_matches(self, config1, config2):
        """
        Check if two configurations match.
        
        Args:
            config1: First config dictionary
            config2: Second config dictionary
            
        Returns:
            bool: True if configs match, False otherwise
        """
        if config1 is None or config2 is None:
            return False
        
        # Keys to compare (ignore others)
        keys_to_compare = ['model', 'temperature', 'num_papers', 'keywords', 'use_abstracts', 'top_abstracts', 'lof_gof', 'knowledge', 'simple', 'proba', 'gencc']
        
        # Boolean keys where None and False should be treated as equivalent
        boolean_keys = {'lof_gof', 'knowledge', 'simple', 'use_abstracts', 'proba', 'gencc'}
        
        for key in keys_to_compare:
            val1 = config1.get(key)
            val2 = config2.get(key)
            
            # Special handling for boolean keys: treat None and False as equivalent
            if key in boolean_keys:
                # Normalize: None, False, 0 -> False; anything else -> its boolean value
                val1 = bool(val1) if val1 is not None else False
                val2 = bool(val2) if val2 is not None else False
                if val1 != val2:
                    return False
            # Special handling for lists (keywords)
            elif isinstance(val1, list) or isinstance(val2, list):
                # Treat None and empty list as equivalent
                val1 = val1 if val1 else []
                val2 = val2 if val2 else []
                if sorted(val1) != sorted(val2):
                    return False
            # Special handling for top_abstracts: None is equivalent to not having it
            elif key == 'top_abstracts':
                if val1 is not None and val2 is not None and val1 != val2:
                    return False
            else:
                if val1 != val2:
                    return False
        
        return True
    
    # =========================================================================
    # Legacy methods for backward compatibility (novel associations, etc.)
    # =========================================================================
    
    def save_mendelian_run(self, gene_name, model, num_papers, keywords, 
                          prompt_template, articles, llm_response_text, parsed_response):
        """
        LEGACY: Save a Mendelian analysis run to JSON file.
        This method is kept for backward compatibility but should be replaced
        with the new run-based system.
        
        Args:
            gene_name: Name of the gene analyzed
            model: LLM model used
            num_papers: Number of papers retrieved
            keywords: Search keywords used
            prompt_template: Prompt template used
            articles: PubMed articles retrieved
            llm_response_text: Raw LLM response
            parsed_response: Parsed response with all results
            
        Returns:
            Filepath where the run was saved, or None if save failed
        """
        try:
            # Get next run number (legacy flat file system)
            run_num = self._get_next_legacy_run_number(gene_name)
            
            filename = f"{gene_name}_run{run_num:02d}.json"
            filepath = os.path.join(self.runs_dir, filename)

            run_data = {
                "gene_symbol": gene_name,
                "run_parameters": {
                    "model": model,
                    "num_papers": num_papers,
                    "keywords": keywords,
                    "prompt": prompt_template
                },
                "pubmed_articles_sent_to_llm": articles,
                "llm_raw_response": llm_response_text,
                "final_output": parsed_response
            }

            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(run_data, f, indent=4)
            
            print(f"INFO: Saved detailed run output to '{filepath}'")
            return filepath

        except Exception as e:
            print(f"ERROR: Could not save run output JSON for {gene_name}. Error: {e}")
            return None
    
    def _get_next_legacy_run_number(self, gene_name):
        """Get the next available run number for a gene (legacy flat file system)."""
        run_num = 1
        while os.path.exists(os.path.join(self.runs_dir, f"{gene_name}_run{run_num:02d}.json")):
            run_num += 1
        return run_num
    
    def load_run(self, gene_name, run_number):
        """
        LEGACY: Load a specific run by gene name and run number.
        
        Args:
            gene_name: Name of the gene
            run_number: Run number to load
            
        Returns:
            Dictionary with run data, or None if not found
        """
        try:
            filename = f"{gene_name}_run{run_number:02d}.json"
            filepath = os.path.join(self.runs_dir, filename)
            
            if not os.path.exists(filepath):
                return None
            
            with open(filepath, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"ERROR: Could not load run {gene_name} run {run_number}: {e}")
            return None
    
    def list_runs(self, gene_name=None):
        """
        LEGACY: List all runs, optionally filtered by gene name.
        
        Args:
            gene_name: Optional gene name to filter by
            
        Returns:
            List of dictionaries with run information
        """
        try:
            if gene_name:
                pattern = os.path.join(self.runs_dir, f"{gene_name}_run*.json")
            else:
                pattern = os.path.join(self.runs_dir, "*_run*.json")
            
            files = glob.glob(pattern)
            
            runs = []
            for filepath in files:
                filename = os.path.basename(filepath)
                # Parse filename to extract gene and run number
                parts = filename.replace('.json', '').split('_run')
                if len(parts) == 2:
                    runs.append({
                        'gene_name': parts[0],
                        'run_number': int(parts[1]),
                        'filepath': filepath,
                        'filename': filename
                    })
            
            # Sort by gene name and run number
            runs.sort(key=lambda x: (x['gene_name'], x['run_number']))
            return runs
            
        except Exception as e:
            print(f"ERROR: Could not list runs: {e}")
            return []
