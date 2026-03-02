import argparse
import time
import asyncio
from typing import List, Tuple, Dict, Optional, Set
import pandas as pd
import requests
import aiohttp
import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

@dataclass
class ProcessingConfig:
    """Configuration for InterPro processing."""
    concurrent_requests: int = 10
    rate_limit_delay: float = 0.1
    timeout: int = 30
    max_retries: int = 3
    retry_delay: float = 2.0
    chunk_size: int = 100

def parse_input_tsv(input_filepath):
    '''Parses the input TSV file and returns a list of Protein IDs and column name of the 2nd column.'''
    df = pd.read_csv(input_filepath, sep='\t')
    if 'UniProt_ID' not in df.columns:
        raise ValueError("Input TSV file must contain a 'UniProt_ID' column.")
    info_name = df.columns[1] if len(df.columns) > 1 else None
    return df['UniProt_ID'].tolist(), info_name, df

async def get_interpro_pfam_info_async(session: aiohttp.ClientSession, protein_id: str, config: ProcessingConfig) -> Optional[List[Dict]]:
    '''Gets Pfam info for a single protein from InterPro API using async requests.'''
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{protein_id}/"
    
    for attempt in range(config.max_retries):
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=config.timeout)) as response:
                if response.status == 200:
                    data = await response.json()
                    return parse_interpro_response(data, protein_id)
                    
                elif response.status == 204:
                    # No content - protein exists but no Pfam domains
                    return []
                    
                elif response.status == 404:
                    # Protein not found
                    print(f"Protein {protein_id} not found in InterPro")
                    return []
                    
                elif response.status == 503:
                    # Service unavailable - retry
                    if attempt < config.max_retries - 1:
                        print(f"Service unavailable for {protein_id}, retrying... (attempt {attempt + 1})")
                        await asyncio.sleep(config.retry_delay * (attempt + 1))  # Exponential backoff
                        continue
                    else:
                        print(f"Failed to retrieve Pfam info for {protein_id}: HTTP {response.status} (max retries reached)")
                        return None
                else:
                    print(f"Failed to retrieve Pfam info for {protein_id}: HTTP {response.status}")
                    return None
                    
        except asyncio.TimeoutError:
            if attempt < config.max_retries - 1:
                print(f"Timeout for {protein_id}, retrying... (attempt {attempt + 1})")
                await asyncio.sleep(config.retry_delay)
            else:
                print(f"Timeout for {protein_id} (max retries reached)")
                return None
        except Exception as e:
            if attempt < config.max_retries - 1:
                print(f"Request error for {protein_id}: {e}, retrying... (attempt {attempt + 1})")
                await asyncio.sleep(config.retry_delay)
            else:
                print(f"Request error for {protein_id}: {e} (max retries reached)")
                return None
    
    return None

def get_interpro_pfam_info(protein_id: str) -> Optional[List[Dict]]:
    '''Synchronous wrapper for backward compatibility.'''
    config = ProcessingConfig()
    return asyncio.run(get_interpro_pfam_info_single(protein_id, config))

async def get_interpro_pfam_info_single(protein_id: str, config: ProcessingConfig) -> Optional[List[Dict]]:
    '''Single protein async request.'''
    connector = aiohttp.TCPConnector(limit=1)
    async with aiohttp.ClientSession(connector=connector) as session:
        return await get_interpro_pfam_info_async(session, protein_id, config)

def parse_interpro_response(data: dict, protein_id: str) -> List[Dict]:
    '''Parse InterPro API response to extract Pfam information.'''
    pfam_info = []
    
    # The response has 'results' array at top level
    results = data.get('results', [])
    
    for result in results:
        # Extract Pfam metadata
        metadata = result.get('metadata', {})
        pfam_id = metadata.get('accession')
        pfam_name = metadata.get('name', '')
        
        if not pfam_id:
            continue
        
        # Look for the specific protein in the proteins array
        proteins = result.get('proteins', [])
        protein_found = False
        best_evalue = None
        
        for protein_entry in proteins:
            # Check if this is our protein (case-insensitive comparison)
            if protein_entry.get('accession', '').lower() == protein_id.lower():
                protein_found = True
                
                # Get E-value from entry_protein_locations
                locations = protein_entry.get('entry_protein_locations', [])
                for location in locations:
                    score = location.get('score')
                    if score is not None:
                        if best_evalue is None or score < best_evalue:
                            best_evalue = score
                break
        
        if protein_found:
            evalue_str = f"{best_evalue:.2e}" if best_evalue is not None else "N/A"
            pfam_info.append({
                'pfam_id': pfam_id,
                'pfam_name': pfam_name,
                'e_value': evalue_str
            })
    
    return pfam_info

async def process_protein_batch(session: aiohttp.ClientSession, protein_ids: List[str], config: ProcessingConfig, start_idx: int = 0) -> Dict[str, Optional[List[Dict]]]:
    '''Process a batch of proteins concurrently.'''
    semaphore = asyncio.Semaphore(config.concurrent_requests)
    
    async def process_single_protein(protein_id: str, idx: int) -> Tuple[str, Optional[List[Dict]]]:
        async with semaphore:
            if idx > 0 and idx % 10 == 0:
                print(f"  Processing batch item {idx + 1}...")
            result = await get_interpro_pfam_info_async(session, protein_id, config)
            await asyncio.sleep(config.rate_limit_delay)  # Rate limiting
            return protein_id, result
    
    tasks = [process_single_protein(protein_id, start_idx + i) for i, protein_id in enumerate(protein_ids)]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    interpro_data = {}
    for result in results:
        if isinstance(result, Exception):
            print(f"Error processing protein: {result}")
            continue
        protein_id, pfam_info = result
        interpro_data[protein_id] = pfam_info
    
    return interpro_data

async def batch_pull_interpro_pfam_async(protein_ids: List[str], config: ProcessingConfig) -> dict:
    '''Pulls Pfam info for proteins using async requests with concurrent processing.'''
    total_proteins = len(protein_ids)
    interpro_data = {}
    
    print(f"Processing {total_proteins} proteins with {config.concurrent_requests} concurrent requests...")
    
    # Create connector with connection pooling
    connector = aiohttp.TCPConnector(
        limit=config.concurrent_requests * 2,
        limit_per_host=config.concurrent_requests
    )
    
    async with aiohttp.ClientSession(connector=connector) as session:
        # Process in chunks to manage memory and provide progress updates
        for i in range(0, total_proteins, config.chunk_size):
            chunk = protein_ids[i:i + config.chunk_size]
            chunk_end = min(i + config.chunk_size, total_proteins)
            
            print(f"Processing chunk {i+1}-{chunk_end}/{total_proteins}...")
            
            chunk_data = await process_protein_batch(session, chunk, config, i)
            interpro_data.update(chunk_data)
            
            print(f"Completed {chunk_end}/{total_proteins} proteins")
    
    return interpro_data

def batch_pull_interpro_pfam(protein_ids: List[str], config: ProcessingConfig = None) -> dict:
    '''Synchronous wrapper for async batch processing.'''
    if config is None:
        config = ProcessingConfig()
    
    return asyncio.run(batch_pull_interpro_pfam_async(protein_ids, config))

def format_interpro_data_for_tsv(interpro_data: dict, original_df: pd.DataFrame) -> pd.DataFrame:
    '''Convert InterPro data to DataFrame format and merge with original data.'''
    
    # Create a copy of the original DataFrame
    result_df = original_df.copy()
    
    # Add Pfam columns if they don't exist
    if 'Pfam_domains' not in result_df.columns:
        result_df['Pfam_domains'] = ''
    if 'Pfam_evalues' not in result_df.columns:
        result_df['Pfam_evalues'] = ''
    
    # Update each protein's Pfam information
    for protein_id, pfam_info in interpro_data.items():
        # Find the row for this protein
        mask = result_df['UniProt_ID'] == protein_id
        
        if pfam_info is None:
            # Failed to retrieve data - keep existing values or set empty
            if mask.any():
                idx = result_df[mask].index[0]
                if pd.isna(result_df.loc[idx, 'Pfam_domains']):
                    result_df.loc[idx, 'Pfam_domains'] = 'FAILED'
                    result_df.loc[idx, 'Pfam_evalues'] = 'FAILED'
                    
        elif len(pfam_info) == 0:
            # No Pfam domains found
            if mask.any():
                idx = result_df[mask].index[0]
                result_df.loc[idx, 'Pfam_domains'] = ''
                result_df.loc[idx, 'Pfam_evalues'] = ''
                
        else:
            # Format Pfam domains and E-values
            domains = []
            evalues = []
            
            for pfam in pfam_info:
                domains.append(pfam['pfam_id'])
                evalues.append(f"{pfam['pfam_id']}:{pfam['e_value']}")
            
            if mask.any():
                idx = result_df[mask].index[0]
                result_df.loc[idx, 'Pfam_domains'] = ';'.join(domains)
                result_df.loc[idx, 'Pfam_evalues'] = ';'.join(evalues)
    
    return result_df

def load_cache(cache_file: Path) -> Dict[str, Optional[List[Dict]]]:
    '''Load cached results from previous runs.'''
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load cache file {cache_file}: {e}")
    return {}

def save_cache(cache_data: dict, cache_file: Path) -> None:
    '''Save results to cache file.'''
    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump(cache_data, f, indent=2)
    except IOError as e:
        print(f"Warning: Could not save cache file {cache_file}: {e}")

def filter_uncached_proteins(protein_ids: List[str], cache_data: dict) -> Tuple[List[str], Dict[str, Optional[List[Dict]]]]:
    '''Filter out proteins that are already cached.'''
    uncached_proteins = [pid for pid in protein_ids if pid not in cache_data]
    cached_proteins = {pid: cache_data[pid] for pid in protein_ids if pid in cache_data}
    return uncached_proteins, cached_proteins

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pull InterPro info for each protein using the InterPro API and save results to a TSV file.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file containing UniProt IDs.')
    parser.add_argument('--output_file', type=str, help='Output file path (if different from input).')
    parser.add_argument('--test_mode', action='store_true', help='If set, processes only a small subset of data for testing.')
    parser.add_argument('--concurrent_requests', type=int, default=10, help='Number of concurrent requests (default: 10)')
    parser.add_argument('--rate_limit', type=float, default=0.1, help='Delay between requests in seconds (default: 0.1)')
    parser.add_argument('--chunk_size', type=int, default=100, help='Process proteins in chunks of this size (default: 100)')
    parser.add_argument('--cache_file', type=str, help='Cache file to store/load results (default: .cache/interpro_cache.json)')
    parser.add_argument('--no_cache', action='store_true', help='Disable caching')
    args = parser.parse_args()

    # Set up configuration
    config = ProcessingConfig(
        concurrent_requests=args.concurrent_requests,
        rate_limit_delay=args.rate_limit,
        chunk_size=args.chunk_size
    )
    
    print(f"Configuration: {config.concurrent_requests} concurrent requests, {config.rate_limit_delay}s rate limit, chunk size {config.chunk_size}")
    
    print(f"Parsing input TSV file: {args.input_file}")
    protein_ids, info_name, original_df = parse_input_tsv(args.input_file)
    print(f"Found {len(protein_ids)} protein IDs in the input file.")

    if args.test_mode:
        print("Test mode enabled: Processing only the first 60 protein IDs.")
        protein_ids = protein_ids[:60]

    # Set up caching
    interpro_data = {}
    if not args.no_cache:
        cache_file = Path(args.cache_file) if args.cache_file else Path('.cache/interpro_cache.json')
        cache_data = load_cache(cache_file)
        uncached_proteins, cached_proteins = filter_uncached_proteins(protein_ids, cache_data)
        interpro_data.update(cached_proteins)
        
        if cached_proteins:
            print(f"Found {len(cached_proteins)} proteins in cache, processing {len(uncached_proteins)} new proteins")
        protein_ids = uncached_proteins
    
    if protein_ids:  # Only process if there are uncached proteins
        print("Pulling Pfam info from InterPro API...")
        start_time = time.time()
        
        new_data = batch_pull_interpro_pfam(protein_ids, config)
        interpro_data.update(new_data)
        
        processing_time = time.time() - start_time
        
        # Update cache
        if not args.no_cache and new_data:
            all_cache_data = load_cache(cache_file) if not args.cache_file or cache_file.exists() else {}
            all_cache_data.update(new_data)
            save_cache(all_cache_data, cache_file)
            print(f"Updated cache with {len(new_data)} new entries")
    else:
        print("All proteins found in cache, skipping API calls")
        processing_time = 0
    
    # Format data and merge with original DataFrame
    final_df = format_interpro_data_for_tsv(interpro_data, original_df)
    
    # Save results
    output_file = args.output_file or args.input_file
    final_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    all_protein_ids = parse_input_tsv(args.input_file)[0][:10 if args.test_mode else len(parse_input_tsv(args.input_file)[0])]
    successful_count = sum(1 for pid in all_protein_ids if interpro_data.get(pid) is not None)
    domains_found_count = sum(1 for pid in all_protein_ids if interpro_data.get(pid) and len(interpro_data[pid]) > 0)
    
    print(f"\nSummary:")
    print(f"Successfully processed: {successful_count}/{len(all_protein_ids)} proteins")
    print(f"Proteins with Pfam domains: {domains_found_count}/{len(all_protein_ids)} proteins")
    print(f"Results saved to: {output_file}")
    if processing_time > 0:
        print(f"Processing completed in {processing_time:.2f} seconds.")
        print(f"Average time per protein: {processing_time/max(len(protein_ids), 1):.3f} seconds")