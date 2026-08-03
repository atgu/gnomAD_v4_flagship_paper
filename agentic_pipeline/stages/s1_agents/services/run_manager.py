"""Service for managing runs with prompt/config comparison."""
from services.run_persistence_service import RunPersistenceService
from services.prompt_collector import PromptCollector


class RunManager:
    """Manages run creation and comparison logic."""
    
    def __init__(self, runs_dir, search_dirs=None):
        """
        Initialize the run manager.
        
        Args:
            runs_dir: Primary directory where runs are stored (for saving)
            search_dirs: Optional list of additional directories to search for cached results
        """
        self.persistence = RunPersistenceService(runs_dir)
        self.prompt_collector = PromptCollector()
        self.search_dirs = search_dirs or []
    
    def find_or_create_run(self, config):
        """
        Find the latest run that matches current prompts/config, or create a new one.
        Searches in primary dir and all search_dirs.
        
        Args:
            config: Dictionary with configuration parameters:
                - model: LLM model
                - temperature: Temperature setting
                - num_papers: Number of papers
                - keywords: List of keywords
                - use_abstracts: Boolean
                - lof_gof: Boolean
        
        Returns:
            tuple: (read_path, write_path) where:
                - read_path: Path to search for cached results (can be agent_runs or website_runs)
                - write_path: Path to write new results (always in primary/website_runs)
        """
        # Collect current prompts
        current_prompts = self.prompt_collector.collect_all_prompts()
        
        read_path = None
        write_path = None
        
        # First, try to find in search_dirs (agent_runs, etc.) for READING
        for search_dir in self.search_dirs:
            result = self._find_matching_run_in_dir(search_dir, current_prompts, config)
            if result:
                print(f"INFO: Found matching run for reading in {search_dir}")
                read_path = result
                break
        
        # Then check primary dir (website_runs)
        result = self._find_matching_run_in_dir(
            self.persistence.runs_dir, current_prompts, config
        )
        if result:
            run_num = result.split('_')[-1]
            print(f"INFO: Found matching run in primary directory: run_{run_num}")
            # If we found in primary, use it for both read and write
            if read_path is None:
                read_path = result
            write_path = result
        else:
            # No match in primary, create new run for writing
            print("INFO: Creating new run in primary directory for writing...")
            write_path = self.persistence.create_new_run(current_prompts, config)
        
        # If we have a read_path from search_dirs but no write_path yet, 
        # we still need a write_path in primary
        if read_path and not write_path:
            # Check if primary has a matching run, otherwise create one
            write_path = self.persistence.create_new_run(current_prompts, config)
        
        # If no read_path found anywhere, use write_path for reading too
        if not read_path:
            read_path = write_path
        
        return read_path, write_path
    
    def _find_matching_run_in_dir(self, runs_dir, current_prompts, config):
        """
        Find a matching run in a specific directory.
        Searches ALL runs in the directory, not just the latest.
        
        Args:
            runs_dir: Directory to search in
            current_prompts: Current prompts dictionary
            config: Current config dictionary
            
        Returns:
            str: Path to matching run, or None if no match
        """
        import os
        import re
        
        # Create a temporary persistence service for this directory
        temp_persistence = RunPersistenceService(runs_dir)
        
        # Get all run directories
        pattern = re.compile(r'^run_(\d+)$')
        run_paths = []
        
        try:
            for item in os.listdir(runs_dir):
                item_path = os.path.join(runs_dir, item)
                if os.path.isdir(item_path) and pattern.match(item):
                    run_paths.append(item_path)
        except Exception:
            return None
        
        if not run_paths:
            return None
        
        # Sort by run number (descending) to check newest first
        run_paths.sort(key=lambda p: int(p.split('_')[-1]), reverse=True)
        
        # Check each run for a match
        for run_path in run_paths:
            run_prompts = temp_persistence.load_run_prompts(run_path)
            run_config = temp_persistence.load_run_config(run_path)
        
            prompts_match = self.prompt_collector.prompts_match(current_prompts, run_prompts)
            config_match = temp_persistence.config_matches(config, run_config)
        
            if prompts_match and config_match:
                run_num = run_path.split('_')[-1]
                print(f"INFO: Found matching run: run_{run_num}")
                return run_path
        
        return None
    
    def load_gene_if_exists(self, run_path, gene_name):
        """
        Try to load an existing gene result from a run.
        
        Args:
            run_path: Path to run directory
            gene_name: Name of gene to load
            
        Returns:
            dict: Gene result if found, None otherwise
        """
        # Create persistence service for this specific run's directory
        # (it might be in agent_runs or website_runs)
        runs_dir = '/'.join(run_path.split('/')[:-1])  # Get parent directory
        temp_persistence = RunPersistenceService(runs_dir)
        
        result = temp_persistence.load_gene_result(run_path, gene_name)
        if result:
            print(f"INFO: Loaded cached result for {gene_name} from {run_path}")
        return result
    
    def save_gene_result(self, run_path, gene_name, result_data):
        """
        Save a gene result to a run.
        
        Args:
            run_path: Path to run directory
            gene_name: Name of gene
            result_data: Analysis results
            
        Returns:
            str: Path to saved file
        """
        # Create persistence service for this specific run's directory
        runs_dir = '/'.join(run_path.split('/')[:-1])  # Get parent directory
        temp_persistence = RunPersistenceService(runs_dir)
        
        return temp_persistence.save_gene_result(run_path, gene_name, result_data)

    def save_gene_raw_llm(self, run_path, gene_name, raw_llm_data):
        """
        Save raw LLM artifacts for a gene to a run (alongside results).
        """
        runs_dir = '/'.join(run_path.split('/')[:-1])  # Get parent directory
        temp_persistence = RunPersistenceService(runs_dir)
        return temp_persistence.save_gene_raw_llm(run_path, gene_name, raw_llm_data)

