from mygene import MyGeneInfo

def get_transcripts_mygene(gene_name):
    """
    Get transcript IDs using MyGene.info service
    
    Args:
        gene_name (str): Gene name (e.g., 'KLHL17')
    
    Returns:
        list: List of Ensembl transcript IDs
    """
    mg = MyGeneInfo()
    results = mg.query(gene_name, 
                      species='human',
                      fields='ensembl.transcript')
    
    # Extract transcript IDs from results
    hits = results.get('hits', [])
    if hits and 'ensembl' in hits[0]:
        transcripts = hits[0]['ensembl']
        if isinstance(transcripts, dict):
            return [transcripts.get('transcript')]
        elif isinstance(transcripts, list):
            return [t.get('transcript') for t in transcripts if 'transcript' in t]
    return []

# Example usage
if __name__ == "__main__":
    gene = "KLHL17"
    transcripts = get_transcripts_mygene(gene)
    print(f"Transcripts for {gene}:")
    for t in transcripts:
        print(t)