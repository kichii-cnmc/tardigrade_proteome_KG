"""
UniProt Protein Information Fetcher

This script fetches protein information from the UniProt API including:
- GO terms (molecular function, cellular component, biological process)
- Pfam domains
- KEGG pathways
- PROSITE annotations
- PDB/AlphaFold structure IDs
- Gene sequences

Supports multiple processing methods: async, batch, and sequential.
"""

import asyncio
import argparse
import time
from typing import Dict, List, Optional, Tuple, Any

import aiohttp
import requests
import pandas as pd

# Constants
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
UNIPROT_SEARCH_URL = f"{UNIPROT_BASE_URL}/search"
UNIPROT_FIELDS = "accession,sequence,go_p,go_c,go_f,go_id,xref_pfam,xref_kegg,xref_prosite,xref_alphafolddb,xref_pdb,xref_geneid"

# Processing parameters
DEFAULT_MAX_CONCURRENT = 10
DEFAULT_BATCH_SIZE = 25
MAX_BATCH_SIZE = 50
DEFAULT_TIMEOUT = 30
BATCH_REQUEST_TIMEOUT = 60
API_DELAY = 0.1
MAX_PDB_STRUCTURES = 5
MAX_RETRIES = 3

# Output column headers
OUTPUT_HEADERS = [
    "UniProt_ID", "NCBI_ID", "Sequence", "GeneID", "GO_mf", "GO_cc", "GO_bp",
    "Pfam_domains", "KEGG_pathways", "PROSITE_annotations", 
    "PDB_structures", "AF_structures"
]

async def fetch_uniprot_info_async(
    session: aiohttp.ClientSession, 
    uniprot_id: str, 
    semaphore: asyncio.Semaphore, 
    max_retries: int = MAX_RETRIES
) -> Optional[Dict[str, Any]]:
    """
    Asynchronously fetches protein information from the UniProt API.
    
    Args:
        session: aiohttp client session
        uniprot_id: UniProt accession ID
        semaphore: Concurrency limiter
        max_retries: Maximum number of retry attempts
        
    Returns:
        Parsed protein information dictionary or None if failed
    """
    params = {"fields": UNIPROT_FIELDS}
    base_url = f"{UNIPROT_BASE_URL}/{uniprot_id}.json"
    
    async with semaphore:  # Limit concurrent requests
        for attempt in range(max_retries):
            try:
                timeout = aiohttp.ClientTimeout(total=DEFAULT_TIMEOUT)
                async with session.get(base_url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return parse_uniprot_data(data, uniprot_id)
                    elif response.status == 429:  # Rate limited
                        wait_time = 2 ** attempt  # Exponential backoff
                        print(f"Rate limited for {uniprot_id}, waiting {wait_time}s...")
                        await asyncio.sleep(wait_time)
                        continue
                    else:
                        if attempt == max_retries - 1:
                            print(f"Failed to fetch data for {uniprot_id}: HTTP {response.status}")
                        continue
            except asyncio.TimeoutError:
                if attempt == max_retries - 1:
                    print(f"Timeout after {DEFAULT_TIMEOUT}s for {uniprot_id}")
                await asyncio.sleep(1)
            except Exception as e:
                if attempt == max_retries - 1:
                    print(f"Unexpected error for {uniprot_id}: {type(e).__name__}: {e}")
                await asyncio.sleep(1)
    
    return None

def fetch_uniprot_batch(uniprot_ids_batch: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Fetches multiple proteins at once using UniProt's batch API.
    
    Args:
        uniprot_ids_batch: List of UniProt IDs to fetch
        
    Returns:
        Dictionary mapping UniProt IDs to their parsed information
    """
    if not uniprot_ids_batch:
        return {}
    
    # Limit batch size to avoid 400 errors
    if len(uniprot_ids_batch) > MAX_BATCH_SIZE:
        uniprot_ids_batch = uniprot_ids_batch[:MAX_BATCH_SIZE]
        print(f"Warning: Batch size limited to {MAX_BATCH_SIZE} entries")
    
    # Build query string for batch request
    ids_str = " OR ".join([f"accession:{uid}" for uid in uniprot_ids_batch])
    params = {
        "query": ids_str,
        "fields": UNIPROT_FIELDS,
        "format": "json",
        "size": len(uniprot_ids_batch)
    }
    
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=BATCH_REQUEST_TIMEOUT)
        
        if response.status_code == 200:
            data = response.json()
            results = {}
            
            for entry in data.get('results', []):
                primary_accession = entry.get('primaryAccession', '')
                if primary_accession:
                    results[primary_accession] = parse_uniprot_data(entry, primary_accession)
                    
            print(f"Successfully fetched {len(results)}/{len(uniprot_ids_batch)} proteins in batch")
            return results
            
        else:
            print(f"Batch request failed: HTTP {response.status_code}")
            if response.status_code == 400:
                print(f"Bad request - query may be too complex. Query length: {len(ids_str)} chars")
            elif response.status_code == 429:
                print("Rate limited - consider reducing batch size or adding delays")
            return {}
            
    except requests.exceptions.Timeout:
        print(f"Batch request timed out after {BATCH_REQUEST_TIMEOUT}s")
        return {}
    except Exception as e:
        print(f"Batch request error: {type(e).__name__}: {e}")
        return {}

def create_empty_protein_info() -> Dict[str, Any]:
    """
    Creates an empty protein information dictionary with all expected fields.
    
    Returns:
        Dictionary with empty sets for all annotation fields and empty string for sequence
    """
    return {
        "GeneID": set(),
        "GO_mf": set(),  # Gene Ontology - Molecular Function
        "GO_cc": set(),  # Gene Ontology - Cellular Component  
        "GO_bp": set(),  # Gene Ontology - Biological Process
        "Pfam_domains": set(),
        "KEGG_pathways": set(),
        "PROSITE_annotations": set(),
        "PDB_structures": set(),
        "AF_structures": set(),  # AlphaFold structures
        "Sequence": ""
    }


def parse_uniprot_data(data: Dict[str, Any], primary_id: str) -> Dict[str, Any]:
    """
    Parses UniProt JSON data to extract relevant information.
    
    Args:
        data: UniProt JSON response data
        primary_id: UniProt accession ID for error reporting
        
    Returns:
        Dictionary containing parsed protein information
    """
    info = create_empty_protein_info()
    
    # Extract cross-references and sequence
    _parse_cross_references(data, info)
    _parse_sequence(data, info)
    
    # Convert sets to lists (except Sequence which remains a string)
    for key, value in info.items():
        if key != "Sequence" and isinstance(value, set):
            info[key] = list(value)
            
    return info


def _parse_cross_references(data: Dict[str, Any], info: Dict[str, Any]) -> None:
    """
    Parse cross-references from UniProt data.
    
    Args:
        data: UniProt JSON response data
        info: Protein information dictionary to populate
    """
    if 'uniProtKBCrossReferences' in data:
        for db_ref in data['uniProtKBCrossReferences']:
            db_type = db_ref.get('database', '')
            ref_id = db_ref.get('id', '')
            
            if not ref_id:
                continue
                
            if db_type == 'GeneID':
                info["GeneID"].add(ref_id)
                
            elif db_type == 'GO':
                _parse_go_term(db_ref, info)
                
            elif db_type == 'Pfam':
                info["Pfam_domains"].add(ref_id)
                
            elif db_type == 'KEGG':
                info["KEGG_pathways"].add(ref_id)
                
            elif db_type == 'PROSITE':
                info["PROSITE_annotations"].add(ref_id)
                
            elif db_type == 'PDB':
                info["PDB_structures"].add(ref_id)
                
            elif db_type == 'AlphaFoldDB':
                info["AF_structures"].add(ref_id)


def _parse_go_term(db_ref: Dict[str, Any], info: Dict[str, Any]) -> None:
    """
    Parse GO term and classify into molecular function, cellular component, or biological process.
    
    Args:
        db_ref: Database reference dictionary containing GO information
        info: Protein information dictionary to populate
    """
    go_id = db_ref.get('id', '')
    properties = db_ref.get('properties', [])
    
    if not properties:
        return
        
    go_description = properties[0].get('value', '')
    
    if not go_description:
        return
        
    # Classify GO term by first character of description
    go_entry = f"{go_id} - {go_description}"
    
    if go_description.startswith('F:'):
        info["GO_mf"].add(go_entry)  # Molecular Function
    elif go_description.startswith('C:'):
        info["GO_cc"].add(go_entry)  # Cellular Component
    elif go_description.startswith('P:'):
        info["GO_bp"].add(go_entry)  # Biological Process


def _parse_sequence(data: Dict[str, Any], info: Dict[str, Any]) -> None:
    """
    Parse protein sequence from UniProt data.
    
    Args:
        data: UniProt JSON response data
        info: Protein information dictionary to populate
    """
    if 'sequence' in data and 'value' in data['sequence']:
        info["Sequence"] = data['sequence']['value']

def fetch_uniprot_info(uniprot_id: str) -> Optional[Dict[str, Any]]:
    """
    Synchronous fallback function for fetching protein information.
    
    Args:
        uniprot_id: UniProt accession ID
        
    Returns:
        Parsed protein information dictionary or None if failed
    """
    params = {"fields": UNIPROT_FIELDS}
    url = f"{UNIPROT_BASE_URL}/{uniprot_id}.json"
    
    try:
        response = requests.get(url, params=params, timeout=DEFAULT_TIMEOUT)
        
        if response.status_code == 200:
            data = response.json()
            return parse_uniprot_data(data, uniprot_id)
        else:
            print(f"Failed to fetch data for {uniprot_id}: HTTP {response.status_code}")
            return None
            
    except requests.exceptions.Timeout:
        print(f"Timeout after {DEFAULT_TIMEOUT}s for {uniprot_id}")
        return None
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {type(e).__name__}: {e}")
        return None


def format_protein_row_dict(uniprot_id: str, ncbi_id: str, info: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Formats protein information into a dictionary for DataFrame creation.
    
    Args:
        uniprot_id: UniProt accession ID
        ncbi_id: NCBI ID
        info: Protein information dictionary or None for failed fetches
        
    Returns:
        Dictionary representing a DataFrame row
    """
    if info is None:
        # Return empty row for failed requests
        return {
            "UniProt_ID": uniprot_id,
            "NCBI_ID": ncbi_id,
            "Sequence": "",
            "GeneID": "",
            "GO_mf": "",
            "GO_cc": "",
            "GO_bp": "",
            "Pfam_domains": "",
            "KEGG_pathways": "",
            "PROSITE_annotations": "",
            "PDB_structures": "",
            "AF_structures": ""
        }
    
    return {
        "UniProt_ID": uniprot_id,
        "NCBI_ID": ncbi_id,
        "Sequence": info.get("Sequence", ""),
        "GeneID": ";".join(info.get("GeneID", [])),
        "GO_mf": ";".join(info.get("GO_mf", [])),
        "GO_cc": ";".join(info.get("GO_cc", [])),
        "GO_bp": ";".join(info.get("GO_bp", [])),
        "Pfam_domains": ";".join(info.get("Pfam_domains", [])),
        "KEGG_pathways": ";".join(info.get("KEGG_pathways", [])),
        "PROSITE_annotations": ";".join(info.get("PROSITE_annotations", [])),
        "PDB_structures": ";".join(info.get("PDB_structures", [])[:MAX_PDB_STRUCTURES]),
        "AF_structures": ";".join(info.get("AF_structures", []))
    }


def save_dataframe_to_tsv(df: pd.DataFrame, output_tsv: str) -> None:
    """
    Save DataFrame to TSV file.
    
    Args:
        df: DataFrame containing protein information
        output_tsv: Output TSV file path
    """
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved {len(df)} protein records to {output_tsv}")

def save_dataframe_to_separate_tsvs(df: pd.DataFrame, output_prefix: str) -> None:
    '''
    Save DataFrame to separate TSV files for each column, includes primary ID column.'''
    for column in df.columns:
        output_tsv = f"{output_prefix}_{column}.tsv"
        df_subset = df[[df.columns[0], column]]  # Include primary ID column
        df_subset.to_csv(output_tsv, sep='\t', index=False)
        print(f"Saved {len(df_subset)} records to {output_tsv}")

async def process_proteins_async(
    uniprot_ids: List[str], 
    ncbi_ids: List[str], 
    max_concurrent: int = DEFAULT_MAX_CONCURRENT
) -> pd.DataFrame:
    """
    Process proteins asynchronously with controlled concurrency.
    
    Args:
        uniprot_ids: List of UniProt IDs to process
        ncbi_ids: List of corresponding NCBI IDs  
        max_concurrent: Maximum number of concurrent requests
        
    Returns:
        DataFrame containing protein information
    """
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async with aiohttp.ClientSession() as session:
        # Create tasks for all proteins
        tasks = [
            fetch_uniprot_info_async(session, uniprot_id, semaphore)
            for uniprot_id in uniprot_ids
        ]
        
        results = []
        
        # Process results as they complete
        print(f"Processing {len(tasks)} proteins...")
        completed = 0
        
        for i, task in enumerate(asyncio.as_completed(tasks)):
            info = await task
            uniprot_id = uniprot_ids[i]
            ncbi_id = ncbi_ids[i] if i < len(ncbi_ids) else ""
            
            row_dict = format_protein_row_dict(uniprot_id, ncbi_id, info)
            results.append(row_dict)
            
            completed += 1
            if completed % 10 == 0:
                print(f"Processed {completed}/{len(tasks)} proteins")
    
    return pd.DataFrame(results)


def process_proteins_batch(
    uniprot_ids: List[str], 
    ncbi_ids: List[str], 
    batch_size: int = DEFAULT_BATCH_SIZE
) -> pd.DataFrame:
    """
    Process proteins in batches using UniProt's batch API.
    Uses smaller batch sizes to avoid API 400 errors.
    
    Args:
        uniprot_ids: List of UniProt IDs to process
        ncbi_ids: List of corresponding NCBI IDs
        batch_size: Number of proteins to process in each batch
        
    Returns:
        DataFrame containing protein information
    """
    results = []
    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
    
    print(f"Processing {len(uniprot_ids)} proteins in {total_batches} batches...")
    
    for i in range(0, len(uniprot_ids), batch_size):
        batch_num = (i // batch_size) + 1
        print(f"Processing batch {batch_num}/{total_batches}")
        
        batch_uniprot_ids = uniprot_ids[i:i+batch_size]
        batch_ncbi_ids = ncbi_ids[i:i+batch_size] if i+batch_size <= len(ncbi_ids) else ncbi_ids[i:] + [""] * (i+batch_size-len(ncbi_ids))
        
        # Fetch batch data
        batch_results = fetch_uniprot_batch(batch_uniprot_ids)
        
        # Process each protein in the batch
        for j, uniprot_id in enumerate(batch_uniprot_ids):
            ncbi_id = batch_ncbi_ids[j] if j < len(batch_ncbi_ids) else ""
            info = batch_results.get(uniprot_id)
            
            # Try individual fetch as fallback if batch failed
            if info is None:
                info = fetch_uniprot_info(uniprot_id)
            
            row_dict = format_protein_row_dict(uniprot_id, ncbi_id, info)
            results.append(row_dict)
        
        # Small delay between batches to be respectful to the API
        time.sleep(API_DELAY)
    
    return pd.DataFrame(results)


def parse_input_tsv(input_tsv: str) -> Tuple[List[str], List[str]]:
    """
    Parses the input TSV file to extract UniProt IDs and NCBI IDs.

    Args:
        input_tsv: Path to the input TSV file

    Returns:
        Tuple of (uniprot_ids, ncbi_ids) lists
    """
    uniprot_ids = []
    ncbi_ids = []
    with open(input_tsv, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if parts[0]:  # Assuming UniProt_ID is the first column
                uniprot_ids.append(parts[0])
            if len(parts) > 2 and parts[2]:  # Assuming NCBI_ID is the third column
                ncbi_ids.append(parts[2])
            else:
                ncbi_ids.append("")
    return uniprot_ids, ncbi_ids


def process_proteins_sequential(
    uniprot_ids: List[str], 
    ncbi_ids: List[str]
) -> pd.DataFrame:
    """
    Process proteins sequentially (original method, kept for compatibility).
    
    Args:
        uniprot_ids: List of UniProt IDs to process
        ncbi_ids: List of corresponding NCBI IDs
        
    Returns:
        DataFrame containing protein information
    """
    results = []
    
    print(f"Processing {len(uniprot_ids)} proteins sequentially...")
    
    for i, (uniprot_id, ncbi_id) in enumerate(zip(uniprot_ids, ncbi_ids)):
        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{len(uniprot_ids)} proteins")
        
        info = fetch_uniprot_info(uniprot_id)
        row_dict = format_protein_row_dict(uniprot_id, ncbi_id, info)
        results.append(row_dict)
    
    return pd.DataFrame(results)


def main() -> None:
    """
    Main function to handle command-line arguments and orchestrate protein processing.
    """
    parser = argparse.ArgumentParser(
        description="Fetch protein information from UniProt API",
        epilog="Supports async, batch, and sequential processing methods for optimal performance."
    )
    parser.add_argument("input_tsv", help="Input TSV file with UniProt IDs")
    parser.add_argument("output_tsv", help="Output TSV file to write protein information")
    parser.add_argument(
        "--method", 
        choices=["async", "batch", "sequential"], 
        default="batch",
        help="Processing method (default: batch)"
    )
    parser.add_argument(
        "--max-concurrent", 
        type=int, 
        default=DEFAULT_MAX_CONCURRENT,
        help=f"Maximum concurrent requests for async method (default: {DEFAULT_MAX_CONCURRENT})"
    )
    parser.add_argument(
        "--batch-size", 
        type=int, 
        default=DEFAULT_BATCH_SIZE,
        help=f"Batch size for batch method (default: {DEFAULT_BATCH_SIZE})"
    )
    args = parser.parse_args()

    # Load UniProt and NCBI IDs from input file
    uniprot_ids, ncbi_ids = parse_input_tsv(args.input_tsv)
    
    print(f"Processing {len(uniprot_ids)} proteins using {args.method} method...")

    # Process proteins using selected method
    start_time = time.time()
    
    if args.method == "async":
        try:
            df = asyncio.run(process_proteins_async(
                uniprot_ids, ncbi_ids, args.max_concurrent
            ))
        except ImportError:
            print("Warning: aiohttp not available. Install with: pip install aiohttp")
            print("Falling back to batch method...")
            df = process_proteins_batch(uniprot_ids, ncbi_ids, args.batch_size)
            
    elif args.method == "batch":
        df = process_proteins_batch(uniprot_ids, ncbi_ids, args.batch_size)
        
    else:  # sequential method
        df = process_proteins_sequential(uniprot_ids, ncbi_ids)
    
    # Save DataFrame to TSV
    # save_dataframe_to_tsv(df, args.output_tsv)
    save_dataframe_to_separate_tsvs(df, args.output_tsv.replace('.tsv', ''))
    
    # Report timing
    elapsed_time = time.time() - start_time
    avg_time_per_protein = elapsed_time / len(uniprot_ids) if uniprot_ids else 0
    
    print(f"Processing completed in {elapsed_time:.2f} seconds")
    print(f"Average time per protein: {avg_time_per_protein:.3f} seconds")


if __name__ == "__main__":
    main()