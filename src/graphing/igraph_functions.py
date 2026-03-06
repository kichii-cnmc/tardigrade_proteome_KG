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


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Predict links in a knowledge graph using R-GCN")
    argparser.add_argument("graphml_file", type=str, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, default="Dsup", help="Comma-separated list of query node IDs (e.g., 'protein1,trait1')")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of predicted top linking proteins to output")
    args = argparser.parse_args()

    # open graph from graphml file with igraph
    g = ig.Graph.Read_GraphML(args.graphml_file)

    # extract triples (source, label, target), edge-weights, and node-types from the graph
    triples = []
    edge_weights = []
    node_types = {}
    for edge in g.es:
        source = g.vs[edge.source]["name"]
        target = g.vs[edge.target]["name"]
        label = edge["edge_type"]
        weight = edge["weight"] if "weight" in edge.attributes() else 1.0
        triples.append((source, label, target))
        edge_weights.append(weight)
    for vertex in g.vs:
        node_id = vertex["name"]
        node_type = vertex["type"] if "type" in vertex.attributes() else "unknown"
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
    results = pipeline(
        training=training,
        testing=testing,
        validation=validation,
        model='ComplEx',  # Use ComplEx for directed relations
        model_kwargs={'embedding_dim': 200},  # Higher dim for better representation
        training_kwargs={'num_epochs': 150, 'batch_size': 512},
        loss='SoftplusLoss',  # Soft loss for weak signals
        loss_kwargs={'reduction': 'mean'},
        device='cpu',
    )
    
    print(f"Finished training ComplEx model with soft loss")

    model = results.model
    print("Finished training the model")

    # # create a weighted loss function that incorporates the edge weights into the training process
    # weights_tensor = torch.tensor(edge_weights, dtype=torch.float)
    # loss = MarginRankingLoss(margin=1.0)
    # training_loop = SLCWATrainingLoop(model=model, triples_factory=tf, loss=loss, sample_weights=weights_tensor)
    # print("Initialized training loop with weighted loss function")

    # # train the model
    # training_loop.train(num_epochs=100, batch_size=256, use_tqdm=True)
    # print("Finished training the model")

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
    print(f"Entities: {list(tf.entity_to_id.keys())[:10]}...")  # Show first 10
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
        
        all_entities = [e for e in tf.entity_to_id.keys() if e not in valid_queries]
        
        for query_idx, query_entity in enumerate(valid_queries):
            print(f"Processing query {query_idx+1}/{len(valid_queries)}: {query_entity}")
            
            # Tail Prediction: query_entity -> relation -> target
            tail_predictions = []
            
            for relation in tf.relation_to_id.keys():
                relation_scores = []
                
                for target_entity in all_entities:
                    try:
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
        
        # Step 5: Generate final ranked list
        final_predictions = []
        for entity, rrf_score in sorted(entity_scores.items(), key=lambda x: x[1], reverse=True):
            # Find best contributing relation
            best_relation = max(entity_relations[entity], key=lambda x: x[2])
            avg_rank = sum(entity_ranks[entity]) / len(entity_ranks[entity])
            
            final_predictions.append({
                'target_entity': entity,
                'rrf_score': rrf_score,
                'best_source': best_relation[0],
                'best_relation': best_relation[1], 
                'best_score': best_relation[2],
                'avg_rank': avg_rank,
                'num_sources': len(entity_ranks[entity]),
                'entity_type': node_types.get(entity, 'unknown')
            })
        
        return final_predictions[:top_k]
    
    def check_novelty(predictions, existing_triples):
        """Check if predictions represent novel links"""
        existing_set = set((s, r, t) for s, r, t in existing_triples)
        
        for pred in predictions:
            pred['is_novel'] = True
            # Check if best source-target combination already exists
            best_triple = (pred['best_source'], pred['best_relation'], pred['target_entity'])
            if best_triple in existing_set:
                pred['is_novel'] = False
        
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
        print(f"     RRF Score: {pred['rrf_score']:.4f} | Avg Rank: {pred['avg_rank']:.1f}")
        print(f"     Best Source: {pred['best_source']} --[{pred['best_relation']}]--> {pred['target_entity']}")
        print(f"     Evidence from {pred['num_sources']} sources | {novelty_flag}")
        print()
    
    # Summary statistics
    novel_count = sum(1 for p in predictions if p['is_novel'])
    protein_count = sum(1 for p in predictions if 'protein' in p['entity_type'].lower())
    
    print(f"=== Summary ===")
    print(f"Novel predictions: {novel_count}/{len(predictions)}")
    print(f"Protein targets: {protein_count}/{len(predictions)}")
    print(f"Multi-source evidence (avg): {sum(p['num_sources'] for p in predictions)/len(predictions):.1f} sources per prediction")
