# call the UniProt API to pull GO terms, Pfam domains, KEGG pathways, PROSITE annotations, PDB/AF structure ID
import requests
import argparse
import time
import asyncio
import aiohttp
from concurrent.futures import ThreadPoolExecutor
import sys
from tqdm import tqdm

async def fetch_uniprot_info_async(session, uniprot_id, semaphore, max_retries=3):
    """
    Asynchronously fetches protein information from the UniProt API.
    """
    params = {
        "fields": "sequence,go_p,go_c,go,go_f,go_id,xref_pfam,xref_kegg,xref_prosite,xref_alphafolddb,xref_pdb,xref_geneid"
    }

    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    async with semaphore:  # Limit concurrent requests
        for attempt in range(max_retries):
            try:
                async with session.get(base_url, params=params, timeout=30) as response:
                    if response.status == 200:
                        data = await response.json()
                        return parse_uniprot_data(data, uniprot_id)
                    elif response.status == 429:  # Rate limited
                        await asyncio.sleep(2 ** attempt)  # Exponential backoff
                        continue
                    else:
                        if attempt == max_retries - 1:
                            print(f"Failed to fetch data for {uniprot_id}: {response.status}")
                        continue
            except asyncio.TimeoutError:
                if attempt == max_retries - 1:
                    print(f"Timeout for {uniprot_id}")
                await asyncio.sleep(1)
            except Exception as e:
                if attempt == max_retries - 1:
                    print(f"Error for {uniprot_id}: {e}")
                await asyncio.sleep(1)
    
    return None

def fetch_uniprot_batch(uniprot_ids_batch):
    """
    Fetches multiple proteins at once using UniProt's batch API.
    """
    if not uniprot_ids_batch:
        return {}
    
    # Limit batch size to avoid 400 errors
    if len(uniprot_ids_batch) > 50:
        uniprot_ids_batch = uniprot_ids_batch[:50]
    
    # Use the correct UniProt API endpoint for batch queries
    ids_str = " OR ".join([f"accession:{uid}" for uid in uniprot_ids_batch])
    params = {
        "query": ids_str,
        "fields": "accession,sequence,go_p,go_c,go_f,go_id,xref_pfam,xref_kegg,xref_prosite,xref_alphafolddb,xref_pdb,xref_geneid",
        "format": "json",
        "size": len(uniprot_ids_batch)
    }
    
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    try:
        response = requests.get(url, params=params, timeout=60)
        if response.status_code == 200:
            data = response.json()
            results = {}
            for entry in data.get('results', []):
                uniprot_id = entry.get('primaryAccession', '')
                if uniprot_id:
                    results[uniprot_id] = parse_uniprot_data(entry, uniprot_id)
            return results
        else:
            print(f"Batch request failed: {response.status_code}")
            if response.status_code == 400:
                print(f"Bad request - possibly too many IDs. Query: {ids_str[:200]}...")
            return {}
    except Exception as e:
        print(f"Batch request error: {e}")
        return {}

def parse_uniprot_data(data, uniprot_id):
    """
    Parses UniProt JSON data to extract relevant information.
    """
    info = {
        "GeneID": set(),
        "GO_mf": set(),
        "GO_cc": set(),
        "GO_bp": set(),
        "Pfam_domains": set(),
        "KEGG_pathways": set(),
        "PROSITE_annotations": set(),
        "PDB_structures": set(),
        "AF_structures": set(),
        "Sequence": ""
    }
    
    # Extract cross-references
    if 'uniProtKBCrossReferences' in data:
        for db_ref in data['uniProtKBCrossReferences']:
            db_type = db_ref.get('database', '')
            if db_type == 'GeneID':
                info["GeneID"].add(db_ref['id'])
            elif db_type == 'GO':
                GO_number = db_ref['id']
                GO_phrase = db_ref.get('properties', [])[0]['value'] if db_ref.get('properties') else ''
                if GO_phrase and GO_phrase[0] == 'F':
                    info["GO_mf"].add(GO_number + " - " + GO_phrase)
                elif GO_phrase and GO_phrase[0] == 'C':
                    info["GO_cc"].add(GO_number + " - " + GO_phrase)
                elif GO_phrase and GO_phrase[0] == 'P':
                    info["GO_bp"].add(GO_number + " - " + GO_phrase)
            elif db_type == 'Pfam':
                info["Pfam_domains"].add(db_ref['id'])
            elif db_type == 'KEGG':
                info["KEGG_pathways"].add(db_ref['id'])
            elif db_type == 'PROSITE':
                info["PROSITE_annotations"].add(db_ref['id'])
            elif db_type == 'PDB':
                info["PDB_structures"].add(db_ref['id'])
            elif db_type == 'AlphaFoldDB':
                info["AF_structures"].add(db_ref['id'])
    
    # Extract sequence
    if 'sequence' in data and 'value' in data['sequence']:
        info["Sequence"] = data['sequence']['value']

    # Convert sets to lists (except Sequence which should remain a string)
    for key in info:
        if key != "Sequence" and isinstance(info[key], set):
            info[key] = list(info[key])
    return info

def fetch_uniprot_info(uniprot_id):
    """
    Synchronous fallback function for fetching protein information.
    """
    params = {
        "fields": "sequence,go_p,go_c,go,go_f,go_id,xref_pfam,xref_kegg,xref_prosite,xref_alphafolddb,xref_pdb,xref_geneid"
    }

    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code == 200:
            data = response.json()
            return parse_uniprot_data(data, uniprot_id)
        else:
            print(f"Failed to fetch data for {uniprot_id}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error for {uniprot_id}: {e}")
        return None

async def process_proteins_async(uniprot_ids, ncbi_ids, output_tsv, max_concurrent=10):
    """
    Process proteins asynchronously with controlled concurrency.
    """
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async with aiohttp.ClientSession() as session:
        tasks = []
        for uniprot_id in uniprot_ids:
            task = fetch_uniprot_info_async(session, uniprot_id, semaphore)
            tasks.append(task)
        
        with open(output_tsv, 'w') as out_f:
            header = ["UniProt_ID", "NCBI_ID", "Sequence", "GeneID", "GO_mf", "GO_cc", "GO_bp", "Pfam_domains", "KEGG_pathways", "PROSITE_annotations", "PDB_structures", "AF_structures"]
            out_f.write('\t'.join(header) + '\n')
            
            # Process results as they complete
            for i, task in enumerate(tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Processing proteins")):
                info = await task
                uniprot_id = uniprot_ids[i]
                ncbi_id = ncbi_ids[i] if i < len(ncbi_ids) else ""
                
                if info is None:
                    # Write empty row for failed requests
                    row = [uniprot_id, ncbi_id, "", "", "", "", "", "", "", "", "", ""]
                else:
                    row = [
                        uniprot_id,
                        ncbi_id,
                        info["Sequence"],
                        ';'.join(info["GeneID"]),
                        ';'.join(info["GO_mf"]),
                        ';'.join(info["GO_cc"]),
                        ';'.join(info["GO_bp"]),
                        ';'.join(info["Pfam_domains"]),
                        ';'.join(info["KEGG_pathways"]),
                        ';'.join(info["PROSITE_annotations"]),
                        ';'.join(info["PDB_structures"][:5]),  # Limit to first 5 PDB structures
                        ';'.join(info["AF_structures"])
                    ]
                out_f.write('\t'.join(row) + '\n')

def process_proteins_batch(uniprot_ids, ncbi_ids, output_tsv, batch_size=25):
    """
    Process proteins in batches using UniProt's batch API.
    Uses smaller batch sizes to avoid API 400 errors.
    """
    with open(output_tsv, 'w') as out_f:
        header = ["UniProt_ID", "NCBI_ID", "Sequence", "GeneID", "GO_mf", "GO_cc", "GO_bp", "Pfam_domains", "KEGG_pathways", "PROSITE_annotations", "PDB_structures", "AF_structures"]
        out_f.write('\t'.join(header) + '\n')
        
        total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
        
        for i in tqdm(range(0, len(uniprot_ids), batch_size), total=total_batches, desc="Processing batches"):
            batch_uniprot_ids = uniprot_ids[i:i+batch_size]
            batch_ncbi_ids = ncbi_ids[i:i+batch_size] if i+batch_size <= len(ncbi_ids) else ncbi_ids[i:] + [""] * (i+batch_size-len(ncbi_ids))
            
            # Fetch batch data
            batch_results = fetch_uniprot_batch(batch_uniprot_ids)
            
            # Write results
            for j, uniprot_id in enumerate(batch_uniprot_ids):
                ncbi_id = batch_ncbi_ids[j] if j < len(batch_ncbi_ids) else ""
                info = batch_results.get(uniprot_id)
                
                if info is None:
                    # Try individual fetch as fallback
                    info = fetch_uniprot_info(uniprot_id)
                    if info is None:
                        row = [uniprot_id, ncbi_id, "", "", "", "", "", "", "", "", "", ""]
                    else:
                        print(f"{info["Sequence"]}")
                        row = [
                            uniprot_id,
                            ncbi_id,
                            info["Sequence"],
                            ';'.join(info["GeneID"]),
                            ';'.join(info["GO_mf"]),
                            ';'.join(info["GO_cc"]),
                            ';'.join(info["GO_bp"]),
                            ';'.join(info["Pfam_domains"]),
                            ';'.join(info["KEGG_pathways"]),
                            ';'.join(info["PROSITE_annotations"]),
                            ';'.join(info["PDB_structures"][:5]),
                            ';'.join(info["AF_structures"])
                        ]
                else:
                    row = [
                        uniprot_id,
                        ncbi_id,
                        info["Sequence"],
                        ';'.join(info["GeneID"]),
                        ';'.join(info["GO_mf"]),
                        ';'.join(info["GO_cc"]),
                        ';'.join(info["GO_bp"]),
                        ';'.join(info["Pfam_domains"]),
                        ';'.join(info["KEGG_pathways"]),
                        ';'.join(info["PROSITE_annotations"]),
                        ';'.join(info["PDB_structures"][:5]),
                        ';'.join(info["AF_structures"])
                    ]
                out_f.write('\t'.join(row) + '\n')
            
            # Small delay between batches to be respectful to the API
            time.sleep(0.1)
def parse_input_tsv(input_tsv):
    """
    Parses the input TSV file to extract UniProt IDs.

    Parameters:
    input_tsv (str): Path to the input TSV file.

    Returns:
    tuple: Lists of UniProt IDs and NCBI IDs.
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch protein information from UniProt.")
    parser.add_argument("input_tsv", help="Input TSV file with UniProt IDs.")
    parser.add_argument("output_tsv", help="Output TSV file to write protein information.")
    parser.add_argument("--method", choices=["async", "batch", "sequential"], default="batch",
                        help="Processing method: async (concurrent requests), batch (batch API), or sequential (original)")
    parser.add_argument("--max-concurrent", type=int, default=10, 
                        help="Maximum concurrent requests for async method (default: 10)")
    parser.add_argument("--batch-size", type=int, default=25,
                        help="Batch size for batch method (default: 25)")
    args = parser.parse_args()

    uniprot_ids, ncbi_ids = parse_input_tsv(args.input_tsv)
    # uniprot_ids = ["P01308"] # test with insulin
    # ncbi_ids = [""] * len(uniprot_ids)
    
    print(f"Processing {len(uniprot_ids)} proteins using {args.method} method...")
    
    start_time = time.time()
    
    if args.method == "async":
        # Use asyncio for concurrent processing
        try:
            asyncio.run(process_proteins_async(uniprot_ids, ncbi_ids, args.output_tsv, args.max_concurrent))
        except ImportError:
            print("aiohttp not available. Install with: pip install aiohttp")
            print("Falling back to batch method...")
            process_proteins_batch(uniprot_ids, ncbi_ids, args.output_tsv, args.batch_size)
    elif args.method == "batch":
        # Use batch API
        process_proteins_batch(uniprot_ids, ncbi_ids, args.output_tsv, args.batch_size)
    else:
        # Original sequential method (kept for compatibility)
        with open(args.output_tsv, 'w') as out_f:
            header = ["UniProt_ID", "NCBI_ID", "GO_mf", "GO_cc", "GO_bp", "Pfam_domains", "KEGG_pathways", "PROSITE_annotations", "PDB_structures", "AF_structures"]
            out_f.write('\t'.join(header) + '\n')
            for uniprot_id, ncbi_id in tqdm(zip(uniprot_ids, ncbi_ids), total=len(uniprot_ids), desc="Processing proteins"):
                print(f"Fetching data for {uniprot_id}...")
                info = fetch_uniprot_info(uniprot_id)
                if info is None:
                    continue
                row = [
                    uniprot_id,
                    ncbi_id,
                    ';'.join(info["GO_mf"]),
                    ';'.join(info["GO_cc"]),
                    ';'.join(info["GO_bp"]),
                    ';'.join(info["Pfam_domains"]),
                    ';'.join(info["KEGG_pathways"]),
                    ';'.join(info["PROSITE_annotations"]),
                    ';'.join(info["PDB_structures"][:5]),  # Limit to first 5 PDB structures
                    ';'.join(info["AF_structures"])
                ]
                out_f.write('\t'.join(row) + '\n')
                # Removed the sleep to make it faster
    
    elapsed_time = time.time() - start_time
    print(f"Processing completed in {elapsed_time:.2f} seconds ({elapsed_time/len(uniprot_ids):.2f} seconds per protein)")