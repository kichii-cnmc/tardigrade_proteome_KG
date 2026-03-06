# script used to make predictions using the KG built by the builder

import random
import torch
import igraph as ig
import numpy as np
import argparse
import pykeen
from pykeen.triples import TriplesFactory
from pykeen.training import SLCWATrainingLoop
from pykeen.losses import MarginRankingLoss
from pykeen.evaluation import RankBasedEvaluator
from pykeen.pipeline import pipeline
from pykeen.losses import SoftplusLoss
from collections import Counter, defaultdict
from typing import List, Dict, Tuple

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

def build_triples_typing_set(g):
    '''Build a set of valid (source_type, relation_type, target_type) triples based on the graph's node and edge attributes'''
    valid_triple_types = set()
    for edge in g.es:
        relation = edge["edge_type"]
        source_type = g.vs[edge.source]["node_type"] if "node_type" in g.vs[edge.source].attributes() else "unknown"
        target_type = g.vs[edge.target]["node_type"] if "node_type" in g.vs[edge.target].attributes() else "unknown"
        valid_triple_types.add((source_type, relation, target_type))
    return valid_triple_types

def check_triple_type_validity(query_triple, valid_triples_set):
    '''Check if a triple's typing exists in the TriplesFactory based on the node_type and edge_type'''
    return query_triple in valid_triples_set

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Predict links in a knowledge graph using R-GCN")
    argparser.add_argument("graphml_file", type=str, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, default="A0A1D1UC20", help="Comma-separated list of query node IDs (e.g., 'protein1,trait1')")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of predicted top linking proteins to output")
    args = argparser.parse_args()

    # open graph from graphml file with igraph
    g = ig.Graph.Read_GraphML(args.graphml_file)

    # extract triples (source, label, target), edge-weights, and node-types from the graph
    triples = []
    edge_weights = []
    node_types = {}
    valid_triple_types = build_triples_typing_set(g)

    for edge in g.es:
        source = g.vs[edge.source]["name"]
        target = g.vs[edge.target]["name"]
        label = edge["edge_type"]
        weight = edge["weight"] if "weight" in edge.attributes() else 1.0
        triples.append((source, label, target))
        edge_weights.append(weight)
    for vertex in g.vs:
        node_id = vertex["name"]
        node_type = vertex["node_type"] if "node_type" in vertex.attributes() else "unknown"
        node_types[node_id] = node_type

    print(f"Extracted {len(triples)} triples from the graph")
    print(f"Extracted {len(edge_weights)} edge weights from the graph")
    print(f"Extracted {len(node_types)} unique nodes from the graph")

    # use TriplesFactory.from_labeled_triples to create a TriplesFactory for pykeen
    triples_array = np.array(triples)
    tf = TriplesFactory.from_labeled_triples(triples_array)
    print(f"Created TriplesFactory with {tf.num_triples} triples, {tf.num_entities} entities, and {tf.num_relations} relations")

    # Create train/validation/test splits
    training, testing, validation = tf.split([0.7, 0.2, 0.1])  # 70/20/10 split

    # Step 2: Knowledge Graph Embedding with ComplEx model and soft loss
    model = pykeen.models.ComplEx(triples_factory=training, embedding_dim=200)
    training_loop = SLCWATrainingLoop(model=model, triples_factory=training)
    
    # Train the model
    training_loop.train(triples_factory=training,num_epochs=150, batch_size=512, use_tqdm=True)
    print("Finished training the model")

    # evaluate the model
    evaluator = RankBasedEvaluator()
    results = evaluator.evaluate(
        model=model,
        mapped_triples=testing.mapped_triples,  # or use tf.testing if you have a test set
        batch_size=256,
        additional_filter_triples=[training.mapped_triples]  # Filter out training triples
    )
    print("Evaluated the model on the training triples (filtered evaluation)")

    print("Evaluation Results:")
    print(f"Hits@1: {results.get_metric('hits_at_1'):.4f}")
    print(f"Hits@3: {results.get_metric('hits_at_3'):.4f}")
    print(f"Hits@10: {results.get_metric('hits_at_10'):.4f}")
    print(f"MRR (Mean Reciprocal Rank): {results.get_metric('mean_reciprocal_rank'):.4f}")
    print(f"MR (Mean Rank): {results.get_metric('mean_rank'):.4f}")

    # Model summary
    print(f"Model: {model}")
    print(f"Number of parameters: {sum(p.numel() for p in model.parameters())}")
    print(f"Number of entities: {model.num_entities}")
    print(f"Number of relations: {model.num_relations}")

    # Check entity and relation mappings
    print(f"Entities: {list(tf.entity_to_id.keys())[:5]}...")  # Show first 10
    print(f"Relations: {list(tf.relation_to_id.keys())}")

    # Step 3: Multi-Source Querying and Step 4: Rank-Based Evidence Aggregation
    def predict_links_with_rrf(model, tf, query_entities: List[str], top_k=10, k_smoothing=60):
        """Implement multi-source querying with Reciprocal Rank Fusion (RRF)"""
        
        # Filter valid query entities
        valid_queries = [q for q in query_entities if q in tf.entity_to_id]
        if not valid_queries:
            print("No valid query entities found in the knowledge graph")
            return []
        
        print(f"Processing {len(valid_queries)} valid query entities: {valid_queries}")
        
        # Step 3a: Parallel Scoring - collect predictions for each query
        entity_scores = defaultdict(float)
        entity_relations = defaultdict(list)
        entity_ranks = defaultdict(list)
        
        # Check entity and relation mappings
        print(f"Entities: {list(tf.entity_to_id.keys())[:5]}...")  # Show first 5
        print(f"Relations: {list(tf.relation_to_id.keys())}")

        # Debug: Check what node types exist
        type_counts = Counter(node_types.values())
        print(f"\nNode types in knowledge graph:")
        for node_type, count in type_counts.most_common():
            print(f"  {node_type}: {count}")
        print()

        # Filter for protein targets only
        all_entities = [e for e in tf.entity_to_id.keys() 
                       if e not in valid_queries and 'protein' in node_types.get(e, '').lower()]
        print(f"Filtering for protein targets: {len(all_entities)} protein entities found")
        
        for query_idx, query_entity in enumerate(valid_queries):
            print(f"Processing query {query_idx+1}/{len(valid_queries)}: {query_entity}")
            
            # Tail Prediction: query_entity -> relation -> target
            tail_predictions = []
            
            for relation in tf.relation_to_id.keys():
                relation_scores = []
                
                for target_entity in all_entities:
                    try:
                        query_type = node_types.get(query_entity, 'unknown')
                        target_type = node_types.get(target_entity, 'unknown')
                        if not check_triple_type_validity((query_type, relation, target_type), valid_triple_types):
                            # print(f"Skipping invalid triple type: ({query_entity}, {relation}, {target_entity})")
                            continue
                        # print(f"Scoring triple: ({query_entity}, {relation}, {target_entity})")
                        # Score the triple (query, relation, target)
                        triple = torch.tensor([[
                            tf.entity_to_id[query_entity],
                            tf.relation_to_id[relation],
                            tf.entity_to_id[target_entity]
                        ]])
                        
                        score = model.score_hrt(triple).item()
                        relation_scores.append((target_entity, relation, score))
                        
                    except KeyError:
                        continue
                
                # Sort by score for this relation
                relation_scores.sort(key=lambda x: x[2], reverse=True)
                tail_predictions.extend(relation_scores)
            
            # Sort all predictions for this query entity
            tail_predictions.sort(key=lambda x: x[2], reverse=True)
            
            # Step 4: Apply RRF - assign ranks and calculate RRF scores
            for rank, (target_entity, relation, score) in enumerate(tail_predictions):
                # Reciprocal Rank Fusion formula
                rrf_score = 1.0 / (k_smoothing + rank + 1)
                entity_scores[target_entity] += rrf_score
                entity_relations[target_entity].append((query_entity, relation, score, rank+1))
                entity_ranks[target_entity].append(rank + 1)
        
        # Step 5: Generate final ranked list with multiple evidence relations
        final_predictions = []
        for entity, rrf_score in sorted(entity_scores.items(), key=lambda x: x[1], reverse=True):
            # Get top 3-5 contributing relations as evidence
            sorted_relations = sorted(entity_relations[entity], key=lambda x: x[2], reverse=True)
            evidence_relations = sorted_relations[:min(5, len(sorted_relations))]
            avg_rank = sum(entity_ranks[entity]) / len(entity_ranks[entity])
            
            final_predictions.append({
                'target_entity': entity,
                'rrf_score': rrf_score,
                'evidence_relations': evidence_relations,  # List of (source, relation, score, rank) tuples
                'best_source': evidence_relations[0][0] if evidence_relations else 'unknown',
                'best_relation': evidence_relations[0][1] if evidence_relations else 'unknown',
                'best_score': evidence_relations[0][2] if evidence_relations else 0.0,
                'avg_rank': avg_rank,
                'num_sources': len(entity_ranks[entity]),
                'entity_type': node_types.get(entity, 'protein')  # Default to protein since we filtered
            })
        
        return final_predictions[:top_k]
    
    def check_novelty(predictions, existing_triples):
        """Check if predictions represent novel links"""
        existing_set = set((s, r, t) for s, r, t in existing_triples)
        
        for pred in predictions:
            pred['is_novel'] = True
            # Check if any evidence relation already exists
            for source, relation, score, rank in pred['evidence_relations']:
                if (source, relation, pred['target_entity']) in existing_set:
                    pred['is_novel'] = False
                    break
        
        return predictions

    # Step 3-5: Execute multi-source querying and RRF aggregation
    query_node_ids = [q.strip() for q in args.query_nodes.split(',')]
    print(f"\n=== Multi-Source Knowledge Graph Prediction ===")
    print(f"Query nodes: {query_node_ids}")
    print(f"Top-K: {args.top_k}")
    
    predictions = predict_links_with_rrf(
        model=model,
        tf=tf,  # Use original tf that has all entities
        query_entities=query_node_ids,
        top_k=args.top_k,
        k_smoothing=60
    )
    
    # Add novelty checking
    predictions = check_novelty(predictions, triples)
    
    # Step 5: Display results with relational context
    print(f"\n=== Top {args.top_k} Predicted Protein Links ===")
    print(f"Using Reciprocal Rank Fusion across {len(query_node_ids)} query nodes\n")
    
    for i, pred in enumerate(predictions):
        novelty_flag = "🆕 NOVEL" if pred['is_novel'] else "♻️ Known"
        print(f"{i+1:2d}. {pred['target_entity']} ({pred['entity_type']})")
        print(f"     RRF Score: {pred['rrf_score']:.4f} | Avg Rank: {pred['avg_rank']:.1f} | {novelty_flag}")
        print(f"     Evidence Relations ({len(pred['evidence_relations'])}):") 
        
        for j, (source, relation, score, rank) in enumerate(pred['evidence_relations']):
            print(f"       {j+1}. {source} --[{relation}]--> {pred['target_entity']} (score: {score:.4f}, rank: {rank})")
        
        print(f"     Multi-source evidence from {pred['num_sources']} query sources")
        print()
    
    # Summary statistics
    if predictions:
        novel_count = sum(1 for p in predictions if p['is_novel'])
        avg_evidence_relations = sum(len(p['evidence_relations']) for p in predictions) / len(predictions)
        avg_sources = sum(p['num_sources'] for p in predictions) / len(predictions)
        
        print(f"=== Summary ===")
        print(f"Protein targets found: {len(predictions)}")
        print(f"Novel predictions: {novel_count}/{len(predictions)} ({novel_count/len(predictions)*100:.1f}%)")
        print(f"Evidence relations per prediction: {avg_evidence_relations:.1f} (avg)")
        print(f"Multi-source evidence: {avg_sources:.1f} query sources per prediction (avg)")
    else:
        print(f"=== No protein predictions found ===")
        print(f"Try expanding query or checking if proteins exist in the knowledge graph")
