"""Service for collecting all prompts used across the application."""
from utils.helpers import load_prompt


class PromptCollector:
    """Collects all prompts used in the analysis pipeline."""
    
    # List of all prompt files used in the application
    # Prompts in subdirectories use path format (e.g., "normal/deep_analysis_disease_agent")
    # Prompts at root level use just the name (e.g., "agreement_analysis_prompt")
    PROMPT_KEYS = [
        # Agreement agent (at root level)
        "agreement_analysis_prompt",
        
        # Deep analysis agents (in normal/ subdirectory)
        "normal/deep_analysis_disease_agent",
        "normal/deep_analysis_penetrance_agent", 
        "normal/deep_analysis_inheritance_agent",
        "normal/deep_analysis_mechanism_agent",
        "normal/deep_analysis_onset_severity_agent",
        "normal/deep_analysis_pathogenicity_agent",
        "normal/deep_analysis_algorithmic_summary_agent",
        
        # Probabilistic agents (in proba/ subdirectory)
        # These output probability distributions instead of point scores
        "proba/deep_analysis_penetrance_agent",
        "proba/deep_analysis_inheritance_agent",
        "proba/deep_analysis_onset_severity_agent",
    ]
    
    @classmethod
    def collect_all_prompts(cls):
        """
        Collect all prompts currently in use.
        
        Returns:
            dict: Dictionary mapping prompt keys to their content
        """
        prompts = {}
        
        for key in cls.PROMPT_KEYS:
            try:
                prompt_content = load_prompt(key)
                if prompt_content:
                    prompts[key] = prompt_content
                else:
                    # If prompt doesn't exist, store None
                    prompts[key] = None
            except Exception as e:
                print(f"WARNING: Could not load prompt '{key}': {e}")
                prompts[key] = None
        
        return prompts
    
    @classmethod
    def prompts_match(cls, prompts1, prompts2):
        """
        Check if two prompt dictionaries are identical.
        
        Args:
            prompts1: First prompt dictionary
            prompts2: Second prompt dictionary
            
        Returns:
            bool: True if all prompts match, False otherwise
        """
        if prompts1 is None or prompts2 is None:
            return False
        
        # Check if all keys match
        keys1 = set(prompts1.keys())
        keys2 = set(prompts2.keys())
        
        if keys1 != keys2:
            return False
        
        # Check if all prompt contents match
        for key in keys1:
            if prompts1[key] != prompts2[key]:
                return False
        
        return True

