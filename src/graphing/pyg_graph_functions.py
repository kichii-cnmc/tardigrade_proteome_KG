import torch
import torch.nn.functional as F
import igraph as ig
import numpy as np
import argparse
from torch_geometric.nn import RGCNConv
from sklearn.metrics import roc_auc_score

def evaluate_link_prediction(model, edge_index, edge_types, num_nodes, device, test_ratio=0.2):
    """
    Evaluate the trained model using standard link prediction metrics.
    """
    model.eval()
    
    # Split edges into train/test
    num_edges = edge_index.size(1)
    num_test = int(num_edges * test_ratio)
    
    # Random permutation for splitting
    perm = torch.randperm(num_edges)
    test_edges = edge_index[:, perm[:num_test]]
    test_edge_types = edge_types[perm[:num_test]]
    
    # Get embeddings
    with torch.no_grad():
        z = model.encode(edge_index, edge_types)
    
    # Evaluation metrics storage
    mrr_scores = []
    hits_at_1 = []
    hits_at_10 = []
    all_scores = []
    all_labels = []
    
    print("Evaluating model...")
    
    for i in range(num_test):
        src, dst = test_edges[0, i].item(), test_edges[1, i].item()
        rel_type = test_edge_types[i].item()
        
        # Generate all possible destination nodes for this source
        all_dsts = torch.arange(num_nodes, device=device)
        src_tensor = torch.full((num_nodes,), src, dtype=torch.long, device=device)
        rel_tensor = torch.full((num_nodes,), rel_type, dtype=torch.long, device=device)
        
        test_edge_batch = torch.stack([src_tensor, all_dsts], dim=0)
        
        # Get scores for all possible destinations
        with torch.no_grad():
            scores = model.decode(z, test_edge_batch, rel_tensor)
            scores = torch.sigmoid(scores)
        
        # Create labels (1 for true destination, 0 for others)
        labels = torch.zeros(num_nodes, device=device)
        labels[dst] = 1
        
        # Remove self-loops for fair evaluation
        if src < len(scores):
            scores[src] = -float('inf')
        
        # Store for AUC calculation
        all_scores.extend(scores.cpu().numpy())
        all_labels.extend(labels.cpu().numpy())
        
        # Calculate ranking metrics
        _, sorted_indices = torch.sort(scores, descending=True)
        rank = (sorted_indices == dst).nonzero(as_tuple=True)[0].item() + 1
        
        # MRR
        mrr_scores.append(1.0 / rank)
        
        # Hits@K
        hits_at_1.append(1.0 if rank <= 1 else 0.0)
        hits_at_10.append(1.0 if rank <= 10 else 0.0)
        
        if i % 100 == 0:
            print(f"  Evaluated {i}/{num_test} test edges")
    
    # Calculate final metrics
    mrr = np.mean(mrr_scores)
    hits_1 = np.mean(hits_at_1)
    hits_10 = np.mean(hits_at_10)
    
    # Calculate AUC
    from sklearn.metrics import roc_auc_score
    all_scores = np.array(all_scores)
    all_labels = np.array(all_labels)
    
    # Remove infinite values for AUC calculation
    valid_mask = np.isfinite(all_scores)
    auc = roc_auc_score(all_labels[valid_mask], all_scores[valid_mask])
    
    return {
        'MRR': mrr,
        'Hits@1': hits_1,
        'Hits@10': hits_10,
        'AUC': auc
    }

def evaluate_link_prediction_fold(model, train_edge_index, train_edge_types, 
                                 test_edge_index, test_edge_types, num_nodes, device):
    """
    Evaluate model on a specific test fold, ensuring no data leakage.
    """
    model.eval()
    
    # Get embeddings using training edges only
    with torch.no_grad():
        z = model.encode(train_edge_index, train_edge_types)
    
    # Evaluation metrics storage
    mrr_scores = [] # mrr = 1/rank, where rank is the position of the true edge in the sorted list of scores
    hits_at_1 = [] # hits@1 = 1 if the true edge is ranked 1st, else 0
    hits_at_10 = [] # hits@10 = 1 if the true edge is ranked in the top 10, else 0
    all_scores = [] # Store scores for AUC calculation
    all_labels = [] # Store labels for AUC calculation (1 for true edge, 0 for others)
    
    num_test = test_edge_index.size(1)
    print(f"Evaluating {num_test} test edges...")
    
    for i in range(num_test):
        src, dst = test_edge_index[0, i].item(), test_edge_index[1, i].item()
        rel_type = test_edge_types[i].item()
        
        # Generate all possible destination nodes for this source
        all_dsts = torch.arange(num_nodes, device=device)
        src_tensor = torch.full((num_nodes,), src, dtype=torch.long, device=device)
        rel_tensor = torch.full((num_nodes,), rel_type, dtype=torch.long, device=device)
        
        test_edge_batch = torch.stack([src_tensor, all_dsts], dim=0)
        
        # Get scores for all possible destinations
        with torch.no_grad():
            scores = model.decode(z, test_edge_batch, rel_tensor)
            scores = torch.sigmoid(scores)
        
        # Create labels (1 for true destination, 0 for others)
        labels = torch.zeros(num_nodes, device=device)
        labels[dst] = 1
        
        # Remove self-loops and training edges for fair evaluation
        if src < len(scores):
            scores[src] = -float('inf')
            
        # Optional: Remove training edges from consideration (filtered evaluation)
        # This prevents the model from being "rewarded" for predicting training edges
        train_mask = torch.zeros(num_nodes, dtype=torch.bool, device=device)
        for j in range(train_edge_index.size(1)):
            if train_edge_index[0, j] == src and train_edge_types[j] == rel_type:
                train_mask[train_edge_index[1, j]] = True
        scores[train_mask] = -float('inf')
        
        # Store for AUC calculation
        valid_mask = scores != -float('inf')
        all_scores.extend(scores[valid_mask].cpu().numpy())
        all_labels.extend(labels[valid_mask].cpu().numpy())
        
        # Calculate ranking metrics
        _, sorted_indices = torch.sort(scores, descending=True)
        rank = (sorted_indices == dst).nonzero(as_tuple=True)[0].item() + 1
        
        # MRR
        mrr_scores.append(1.0 / rank)
        
        # Hits@K
        hits_at_1.append(1.0 if rank <= 1 else 0.0)
        hits_at_10.append(1.0 if rank <= 10 else 0.0)
        
        if i % 50 == 0:
            print(f"  Evaluated {i}/{num_test} test edges")
    
    # Calculate final metrics
    mrr = np.mean(mrr_scores)
    hits_1 = np.mean(hits_at_1)
    hits_10 = np.mean(hits_at_10)
    
    # Calculate AUC
    if len(all_scores) > 0:
        auc = roc_auc_score(all_labels, all_scores)
    else:
        auc = 0.0
    
    return {
        'MRR': mrr,
        'Hits@1': hits_1,
        'Hits@10': hits_10,
        'AUC': auc
    }

def evaluate_ranking_quality(model, edge_index, edge_types, query_nodes, target_nodes, num_relations, device):
    """
    Evaluate ranking quality for specific query-target pairs.
    """
    model.eval()
    
    with torch.no_grad():
        z = model.encode(edge_index, edge_types)
    
    ranking_metrics = {}
    
    for query_idx in query_nodes:
        query_scores = {}
        
        for rel_idx in range(num_relations):
            # Create pairs: (Query Node, All Target Nodes)
            q_tensor = torch.full((len(target_nodes),), query_idx, dtype=torch.long, device=device)
            target_tensor = torch.tensor(target_nodes, dtype=torch.long, device=device)
            test_edges = torch.stack([q_tensor, target_tensor], dim=0)
            
            # Use specific relation
            r_tensor = torch.full((len(target_nodes),), rel_idx, dtype=torch.long, device=device)
            
            scores = torch.sigmoid(model.decode(z, test_edges, r_tensor))
            query_scores[rel_idx] = scores.cpu().numpy()
        
        ranking_metrics[query_idx] = query_scores
    
    return ranking_metrics

class RGCNLinkPrediction(torch.nn.Module):
    def __init__(self, num_nodes, num_relations, hidden_dim):
        super().__init__()
        # Node embeddings (learnable since we may not have node features)
        self.node_emb = torch.nn.Embedding(num_nodes, hidden_dim)
        
        # R-GCN layers
        self.conv1 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        self.conv2 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        
        # Relation embeddings for DistMult decoder
        self.rel_emb = torch.nn.Embedding(num_relations, hidden_dim)

    def encode(self, edge_index, edge_type):
        x = self.node_emb.weight
        x = self.conv1(x, edge_index, edge_type)
        x = F.relu(x)
        x = F.dropout(x, p=0.2, training=self.training)
        x = self.conv2(x, edge_index, edge_type)
        return x

    def decode(self, z, edge_index, edge_type):
        # DistMult scoring function: sum(z_src * r * z_dst)
        src, dst = edge_index
        r = self.rel_emb(edge_type)
        return (z[src] * r * z[dst]).sum(dim=-1)


def mapping_types_to_indices(g):
    node_type_mapping = {ntype: idx for idx, ntype in enumerate(set(g.vs["node_type"]))} if "node_type" in g.vs.attributes() else None
    edge_type_mapping = {etype: idx for idx, etype in enumerate(set(g.es["edge_type"]))}
    return node_type_mapping, edge_type_mapping

def map_valid_triple_typings(g):
    '''Create a mapping of valid node type - edge type - node type combinations for filtering predictions.'''
    valid_triple_mapping = {}
    for edge in g.es:
        src_type = g.vs[edge.source]["node_type"] if "node_type" in g.vs.attributes() else None
        dst_type = g.vs[edge.target]["node_type"] if "node_type" in g.vs.attributes() else None
        edge_type = edge["edge_type"]
        
        if src_type and dst_type:
            if edge_type not in valid_triple_mapping:
                valid_triple_mapping[edge_type] = set()
            valid_triple_mapping[edge_type].add((src_type, dst_type))
    return valid_triple_mapping

def build_alias_GO_alias_mappings(g):
    '''Build mappings for GO term IDs to their full node names for inputting queries as GO IDs.'''
    go_id_to_name = {}
    for v in g.vs:
        if v["node_type"] == "Biological_Process" or v["node_type"] == "Molecular_Function" or v["node_type"] == "Cellular_Component":
            go_id = v["name"][:10]  # Assuming the 'name' attribute contains the GO ID
            go_id_to_name[go_id] = v["name"]  # Map GO ID to full node name (could be the same)
    print(f"Built GO ID to name mapping for {len(go_id_to_name)} GO terms")
    return go_id_to_name


def five_fold_cross_validation_training(edge_index, edge_types, num_nodes, num_relations, device, hidden_dim=64, epochs=100):
    """
    Perform 5-fold cross validation on edge prediction, outputting average metrics, fold results, and best model.
    """
    from sklearn.model_selection import KFold
    import torch
    
    # Convert to numpy for sklearn KFold
    num_edges = edge_index.size(1)
    edge_indices = np.arange(num_edges)
    
    # Initialize KFold
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    
    fold_results = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(edge_indices)):
        print(f"\n=== FOLD {fold_idx + 1}/5 ===")
        
        # Split edges into train/test for this fold
        train_edges = edge_index[:, train_idx]
        train_edge_types = edge_types[train_idx]
        test_edges = edge_index[:, test_idx]
        test_edge_types = edge_types[test_idx]
        
        print(f"Train edges: {len(train_idx)}, Test edges: {len(test_idx)}")
        
        # Initialize fresh model for this fold
        model = RGCNLinkPrediction(num_nodes, num_relations, hidden_dim).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
        criterion = torch.nn.BCEWithLogitsLoss()
        
        # Train on training edges only
        print("Training model...")
        for epoch in range(1, epochs + 1):
            model.train()
            optimizer.zero_grad()
            
            # Get node embeddings using ONLY training edges
            z = model.encode(train_edges.to(device), train_edge_types.to(device))
            
            # Positive samples (training edges)
            pos_out = model.decode(z, train_edges.to(device), train_edge_types.to(device))
            pos_loss = criterion(pos_out, torch.ones_like(pos_out))
            
            # Negative sampling - corrupt destination nodes
            neg_dst = torch.randint(0, num_nodes, (train_edges.size(1),), device=device)
            neg_edge_index = torch.stack([train_edges[0].to(device), neg_dst], dim=0)
            neg_out = model.decode(z, neg_edge_index, train_edge_types.to(device))
            neg_loss = criterion(neg_out, torch.zeros_like(neg_out))
            
            # Total loss
            loss = pos_loss + neg_loss
            loss.backward()
            optimizer.step()
            
            if epoch % 20 == 0 or epoch == 1:
                print(f"  Epoch {epoch:03d}/{epochs}, Loss: {loss.item():.4f}")
        
        # Evaluate on test edges
        print("Evaluating fold...")
        fold_metrics = evaluate_link_prediction_fold(
            model, train_edges.to(device), train_edge_types.to(device),
            test_edges.to(device), test_edge_types.to(device), 
            num_nodes, device
        )
        
        fold_results.append(fold_metrics)
        print(f"Fold {fold_idx + 1} Results:")
        for metric, value in fold_metrics.items():
            print(f"  {metric}: {value:.4f}")
    
    # Average results across folds
    avg_results = {}
    for metric in fold_results[0].keys():
        avg_results[metric] = np.mean([fold[metric] for fold in fold_results])
        avg_results[f"{metric}_std"] = np.std([fold[metric] for fold in fold_results])
    
    return avg_results, fold_results, fold_results[np.argmax([fold['MRR'] for fold in fold_results])]['model'] if 'model' in fold_results[0] else None

def predict_edges_for_query_nodes(model, g, query_node_indices, edge_index, edge_types, node_type_map, edge_type_map, valid_triple_typings, device, top_k=10):
    """
    Predict edges for specific query nodes, filtering for non-existing edges and valid node-edge-node type combinations.
    """
    all_predictions = {}
    for query_idx in query_node_indices:
        query_name = g.vs[query_idx]["name"] if has_names else f"node_{query_idx}"
        print(f"\nPredicting edges for query node: {query_name}")
        predictions_for_query = []
        # For each relation type, predict edges to all other nodes
        for rel_type in range(num_relations):
            # Create candidate edges: query_node -> all_other_nodes
            target_nodes = list(range(num_nodes))

            # Remove self-loops
            if query_idx in target_nodes:
                target_nodes.remove(query_idx)
            
            if target_nodes:
                # Create edge tensors
                query_tensor = torch.full((len(target_nodes),), query_idx, dtype=torch.long, device=device)
                target_tensor = torch.tensor(target_nodes, dtype=torch.long, device=device)
                candidate_edges = torch.stack([query_tensor, target_tensor], dim=0)
                rel_tensor = torch.full((len(target_nodes),), rel_type, dtype=torch.long, device=device)
                
                # Get prediction scores
                with torch.no_grad():
                    scores = torch.sigmoid(model.decode(z, candidate_edges, rel_tensor))
                
                # Filter out existing edges and collect predictions
                for i, target_idx in enumerate(target_nodes):
                    edge_tuple = (query_idx, target_idx, rel_type)
                    
                    # Skip if edge already exists
                    if edge_tuple in existing_edges:
                        continue
                    
                    # Only allow certain node types to connect via certain edge types
                    # get node types for query and target
                    query_node_type = g.vs[query_idx]["node_type"] if "node_type" in g.vs.attributes() else None
                    target_node_type = g.vs[target_idx]["node_type"] if "node_type" in g.vs.attributes() else None
                    if rel_type in valid_triple_typings:
                        if (query_node_type, target_node_type) not in valid_triple_typings[rel_type]:
                            continue
                    predictions_for_query.append({
                        'target_node': g.vs[target_idx]["name"] if has_names else f"node_{target_idx}",
                        'edge_type': rel_type,
                        'score': scores[i].item()
                    })
        
        # Sort by score and keep top-k
        predictions_for_query.sort(key=lambda x: x['score'], reverse=True)
        top_predictions = predictions_for_query[:args.top_k]
        
        all_predictions[query_name] = top_predictions

    # Print predictions
    for query_name, predictions in all_predictions.items():
        print(f"\nTop-{args.top_k} predictions for {query_name}:")
        for pred in predictions:
            print(f"  Target Node: {pred['target_node']}, Edge Type: {pred['edge_type']}, Score: {pred['score']:.4f}")

    return all_predictions

def rank_proximal_nodes(g, query_node_indices, node_type='Protein', top_k=10):
    """Identify top-k proximal nodes to multiple query nodes using multi-source RWR."""

    # Convert to undirected if needed
    if not g.is_directed():
        print("Graph is already undirected.")
    else:
        print("Converting graph to undirected for proximity ranking...")
        g = g.as_undirected(combine_edges=None)

    # Filter nodes by type
    if "node_type" in g.vs.attributes():
        target_node_indices = [i for i in range(g.vcount()) if g.vs[i]["node_type"] == node_type]
        print(f"Found {len(target_node_indices)} nodes of type '{node_type}'")
    else:
        target_node_indices = list(range(g.vcount()))
        print(f"No node type filtering available, using all nodes")

    if not target_node_indices:
        return {}

    restart_prob = 0.15
    max_iterations = 1000
    tolerance = 1e-6

    # Adjacency
    adj_matrix = np.array(g.get_adjacency().data, dtype=float)
    n_nodes = adj_matrix.shape[0]

    # Transition matrix with self-loop fix
    row_sums = adj_matrix.sum(axis=1)
    transition_matrix = np.zeros_like(adj_matrix)

    for i in range(n_nodes):
        if row_sums[i] > 0:
            transition_matrix[i] = adj_matrix[i] / row_sums[i]
        else:
            transition_matrix[i, i] = 1.0

    # =========================
    # 🔥 MULTI-SOURCE RESTART
    # =========================
    restart_vector = np.zeros(n_nodes)
    valid_queries = [q for q in query_node_indices if row_sums[q] > 0]

    if not valid_queries:
        print("No valid (non-isolated) query nodes.")
        return {}

    for q in valid_queries:
        restart_vector[q] = 1.0 / len(valid_queries)

    prob_vector = restart_vector.copy()

    print(f"\nRunning multi-source RWR on {len(valid_queries)} query nodes")

    # RWR iteration
    for iteration in range(max_iterations):
        prev_prob = prob_vector.copy()

        prob_vector = (
            (1 - restart_prob) * transition_matrix.T @ prob_vector
            + restart_prob * restart_vector
        )

        if np.linalg.norm(prob_vector - prev_prob, ord=1) < tolerance:
            print(f"Converged in {iteration + 1} iterations")
            break
    else:
        print("Warning: did not converge")

    # Debug
    print("Total mass:", prob_vector.sum())
    print("Max value:", prob_vector.max())
    print("Nonzero count:", np.count_nonzero(prob_vector))

    # =========================
    # OUTPUT (single ranking)
    # =========================
    valid_targets = [
        {
            "node_name": g.vs[i]["name"] if "name" in g.vs.attributes() else f"node_{i}",
            "node_idx": i,
            "proximity_score": prob_vector[i],
        }
        for i in target_node_indices
        if i not in query_node_indices
    ]

    top_targets = sorted(valid_targets, key=lambda x: x["proximity_score"], reverse=True)[:top_k]

    print(f"\nTop {len(top_targets)} proximal nodes (combined):")
    for i, t in enumerate(top_targets, 1):
        print(f"{i}. {t['node_name']} ({t['proximity_score']:.6f})")

    return {"combined_query": top_targets}

def interpret_proximity(g, query_node_indices, proximal_nodes, max_neighbors=5):
    """
    Explain proximity to a SET of query nodes (multi-source RWR).
    """

    results = []

    has_edge_type = "edge_type" in g.es.attributes()
    has_name = "name" in g.vs.attributes()

    query_names = [
        g.vs[i]["name"] if has_name else f"node_{i}"
        for i in query_node_indices
    ]

    print(f"\nInterpreting proximity for query set: {', '.join(query_names)}")

    # Precompute neighbors for all queries
    query_neighbors_map = {
        q: set(g.neighbors(q)) for q in query_node_indices
    }

    # Union of all query neighbors
    all_query_neighbors = set().union(*query_neighbors_map.values())

    for target in proximal_nodes.get("combined_query", []):
        target_idx = target["node_idx"]
        target_name = target["node_name"]

        explanation = {
            "target": target_name,
            "connections": []
        }

        target_neighbors = set(g.neighbors(target_idx))

        # =========================
        # 1. DIRECT CONNECTIONS TO QUERY SET
        # =========================
        direct_connections = []
        edge_types_map = {}

        for q in query_node_indices:
            if g.are_adjacent(q, target_idx):
                q_name = g.vs[q]["name"] if has_name else f"node_{q}"
                direct_connections.append(q_name)

                edge_ids = g.es.select(_between=([q], [target_idx]))
                types = []

                for e in edge_ids:
                    if has_edge_type:
                        types.append(e["edge_type"])
                    else:
                        types.append("direct interaction")

                edge_types_map[q_name] = list(set(types))

        if direct_connections:
            if len(direct_connections) == 1:
                q = direct_connections[0]

                explanation["connections"].append(
                    f"Direct connection to node {q} via: {', '.join(edge_types_map[q])}"
                )
            else:
                explanation["connections"].append(
                    f"Directly connected to MULTIPLE query proteins: {', '.join(direct_connections)}"
                )

                for q in direct_connections[:max_neighbors]:
                    explanation["connections"].append(
                        f"  ↳ {q} via: {', '.join(edge_types_map[q])}"
                    )

        # =========================
        # 2. SHARED NEIGHBORS WITH QUERY SET
        # =========================
        shared = target_neighbors & all_query_neighbors

        if shared:
            shared_names = [
                g.vs[i]["name"] if has_name else f"node_{i}"
                for i in list(shared)[:max_neighbors]
            ]

            explanation["connections"].append(
                f"Shares functional neighbors with query set: {', '.join(shared_names)}"
            )

        # =========================
        # 3. BRIDGING PATHS (TARGET CONNECTS QUERIES)
        # =========================
        bridging_queries = []

        for q in query_node_indices:
            q_neighbors = query_neighbors_map[q]
            if target_idx in q_neighbors:
                bridging_queries.append(
                    g.vs[q]["name"] if has_name else f"node_{q}"
                )

        if len(bridging_queries) >= 2:
            explanation["connections"].append(
                f"Acts as a bridge between query proteins: {', '.join(bridging_queries)}"
            )

        # =========================
        # 4. FALLBACK
        # =========================
        if not explanation["connections"]:
            explanation["connections"].append(
                "No strong local connections; proximity likely arises from global network structure"
            )

        results.append(explanation)

        # =========================
        # PRINT
        # =========================
        print(f"\n  Target: {target_name}")
        for c in explanation["connections"]:
            print(f"    - {c}")

    return {"combined_query": results}

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Run R-GCN link prediction on a graphml file")
    argparser.add_argument("graphml_file", type=str, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, required=True, help="List of query node names separated by commas")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of top predictions to output")
    argparser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    args = argparser.parse_args()

    # load graph as igraph object
    print("Loading graph...")
    g = ig.Graph.Read_GraphML(args.graphml_file)

    # initialize tensor data for PyG (edge_index, edge_types, node_types)
    edge_index = torch.tensor(g.get_edgelist(), dtype=torch.long).t().contiguous()
    node_type_map, edge_type_map = mapping_types_to_indices(g) 
    edge_types = torch.tensor([edge_type_map[etype] for etype in g.es["edge_type"]], dtype=torch.long)
    node_types = torch.tensor([node_type_map[ntype] for ntype in g.vs["node_type"]], dtype=torch.long) if node_type_map else None
    valid_triple_typings = map_valid_triple_typings(g)

    # evaluate graph for basic checks (number of nodes, edges, types)
    print(f"Edge Types: {edge_type_map}")
    num_nodes = g.vcount()
    num_relations = int(edge_types.max().item() + 1)
        # check if graph has node names for mapping
    has_names = "name" in g.vs.attributes()
    print(f"Graph loaded: {num_nodes} nodes, {g.ecount()} edges, {num_relations} relations.")

    # node type check
    from collections import Counter
    node_types = g.vs["node_type"] if "node_type" in g.vs.attributes() else []
    print(Counter(node_types))

    # Parse query node names and map to indices
    go_id_to_name = build_alias_GO_alias_mappings(g)
    query_node_names = [name.strip() for name in args.query_nodes.split(',')]
    for name in query_node_names:
        if name in go_id_to_name:
            print(f"Mapping GO ID '{name}' to node name '{go_id_to_name[name]}'")
            query_node_names[query_node_names.index(name)] = go_id_to_name[name]
    print(f"Query node names: {query_node_names}")
    query_node_indices = []

    # initialize model, optimizer, and loss function
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = RGCNLinkPrediction(num_nodes, num_relations, hidden_dim=64).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    criterion = torch.nn.BCEWithLogitsLoss()
    edge_index, edge_types = edge_index.to(device), edge_types.to(device)

    # five-fold cross validation training for link prediction model
    avg_results, fold_results, best_model = five_fold_cross_validation_training(edge_index, edge_types, num_nodes, num_relations, device, hidden_dim=64, epochs=args.epochs)

    # output evaluation metrics (MRR, Hits@1, Hits@10, AUC) & select best model
    print("\nFive-Fold Cross Validation Results:")
    for metric, value in avg_results.items():
        print(f"  {metric}: {value:.3f}")
    if best_model:
        print("Best model selected based on highest MRR.")
        model = best_model

    # use model to predict edges into graph for query nodes, filter for non-existing edges and correct node-edge-node types
    print("\nPredicting edges for query nodes...")

    if has_names: # checks if graph has node names
        name_to_index = {g.vs[i]["name"]: i for i in range(num_nodes)}
        for name in query_node_names:
            if name in name_to_index:
                query_node_indices.append(name_to_index[name])
                print(f"Found query node '{name}' at index {name_to_index[name]}")
            else:
                print(f"Warning: Query node '{name}' not found in graph")
    else:
        # If no names, treat query_nodes as indices
        try:
            query_node_indices = [int(name) for name in query_node_names]
        except ValueError:
            print("Error: Graph has no node names, please provide node indices")
            exit(1)
    if not query_node_indices:
        print("No valid query nodes found")
        exit(1)

    # Get existing edges for filtering
    existing_edges = set()
    for i in range(edge_index.size(1)):
        src, dst, etype = edge_index[0, i].item(), edge_index[1, i].item(), edge_types[i].item()
        existing_edges.add((src, dst, etype))

    # Set model to evaluation mode and get embeddings
    model.eval()
    with torch.no_grad():
        z = model.encode(edge_index, edge_types)

    # Store predictions for each query node
    edge_expansion_k = g.ecount() * 0.02 # add up to 2% of new edges
    all_predictions = predict_edges_for_query_nodes(model, g, query_node_indices, edge_index, edge_types, node_type_map, edge_type_map, valid_triple_typings, device, top_k=10)
    
    # add prediction into graph as new edges
    for query_name, predictions in all_predictions.items():
        query_idx = name_to_index[query_name] if has_names else int(query_name.split('_')[1])
        for pred in predictions:
            target_name = pred['target_node']
            target_idx = name_to_index[target_name] if has_names else int(target_name.split('_')[1])
            edge_type = pred['edge_type']
            score = pred['score']
            # Add edge to graph with predicted edge type and score as weight
            g.add_edge(query_idx, target_idx, edge_type=edge_type, weight=score)

    # based on new graph w predicted edges, identify highest proximity/correlation nodes to query nodes
    proximal_nodes = rank_proximal_nodes(g, query_node_indices, node_type='Protein', top_k=10)

    # output interpretability analysis results (e.g., shared neighbors, path lengths) for top predicted nodes to provide insights into why those nodes are proximal to query nodes
    interpret_proximity(g, query_node_indices, proximal_nodes)