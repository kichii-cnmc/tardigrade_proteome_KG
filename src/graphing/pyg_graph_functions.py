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

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Run R-GCN link prediction on a graphml file")
    argparser.add_argument("--graphml_file", type=str, required=True, help="Path to the input graphml file")
    argparser.add_argument("--query_nodes", type=str, required=True, help="List of query node names separated by commas")
    argparser.add_argument("--top_k", type=int, default=10, help="Number of top predictions to output")
    argparser.add_argument("--target_node_type", type=str, default="Protein", help="Node type to rank (e.g., protein, trait)")
    argparser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    args = argparser.parse_args()

    print("Loading graph...")
    g = ig.Graph.Read_GraphML(args.graphml_file)
    
    edge_index = torch.tensor(g.get_edgelist(), dtype=torch.long).t().contiguous()
    node_type_map, edge_type_map = mapping_types_to_indices(g)
    edge_types = torch.tensor([edge_type_map[etype] for etype in g.es["edge_type"]], dtype=torch.long)
    node_types = torch.tensor([node_type_map[ntype] for ntype in g.vs["node_type"]], dtype=torch.long) if node_type_map else None

    print(f"Edge Types: {edge_type_map}")
    num_nodes = g.vcount()
    num_relations = int(edge_types.max().item() + 1)
    
    # check if graph has node names for mapping
    has_names = "name" in g.vs.attributes()

    print(f"Graph loaded: {num_nodes} nodes, {g.ecount()} edges, {num_relations} relations.")

    # initialize model, optimizer, and loss
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = RGCNLinkPrediction(num_nodes, num_relations, hidden_dim=64).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    criterion = torch.nn.BCEWithLogitsLoss()

    edge_index, edge_types = edge_index.to(device), edge_types.to(device)

    # train model
    print("Training R-GCN model...")
    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad()
        
        # Get node embeddings
        z = model.encode(edge_index, edge_types)
        
        # Positive samples (actual edges)
        pos_out = model.decode(z, edge_index, edge_types)
        pos_loss = criterion(pos_out, torch.ones_like(pos_out))
        
        # Negative samples (corrupt destination nodes randomly)
        neg_dst = torch.randint(0, num_nodes, (edge_index.size(1),), device=device)
        neg_edge_index = torch.stack([edge_index[0], neg_dst], dim=0)
        neg_out = model.decode(z, neg_edge_index, edge_types)
        neg_loss = criterion(neg_out, torch.zeros_like(neg_out))
        
        # Total loss
        loss = pos_loss + neg_loss
        loss.backward()
        optimizer.step()
        
        if epoch % 20 == 0 or epoch == 1:
            print(f"  Epoch {epoch:03d}/{args.epochs}, Loss: {loss.item():.4f}")

    # evaluate model
    # Evaluate model performance
    print("\nEvaluating model performance...")
    eval_metrics = evaluate_link_prediction(model, edge_index, edge_types, num_nodes, device)

    print(f"\nEvaluation Results:")
    print(f"Mean Reciprocal Rank (MRR): {eval_metrics['MRR']:.4f}")
    print(f"Hits@1: {eval_metrics['Hits@1']:.4f}")
    print(f"Hits@10: {eval_metrics['Hits@10']:.4f}")
    print(f"AUC: {eval_metrics['AUC']:.4f}")

    # 3. Filtered Inference
    model.eval()
    with torch.no_grad():
        z = model.encode(edge_index, edge_types)

    # Identify which node indices match the requested target_node_type
    target_type_id = node_type_map.get(args.target_node_type)
    if target_type_id is None:
        print(f"Error: Node type '{args.target_node_type}' not found in graph labels:")
        exit()

    # Find all indices where node_type matches the filter
    valid_target_indices = (node_types == target_type_id).nonzero(as_tuple=True)[0].to(device)

    query_names = [n.strip() for n in args.query_nodes.split(",")]
    for q_name in query_names:
        try:
            q_idx = g.vs.find(name=q_name).index
        except:
            print(f"Node {q_name} not found."); continue

        print(f"\nPredictions for Correlation with {q_name} (Limited to type: {args.target_node_type})")
        results = []
        for r_idx in range(num_relations):
            # Create pairs: (Query Node, All Valid Target Nodes)
            q_tensor = torch.full((len(valid_target_indices),), q_idx, dtype=torch.long, device=device)
            test_edges = torch.stack([q_tensor, valid_target_indices], dim=0)
            
            # Use relation r_idx for all these pairs
            r_tensor = torch.full((len(valid_target_indices),), r_idx, dtype=torch.long, device=device)
            
            probs = torch.sigmoid(model.decode(z, test_edges, r_tensor))
            
            for i, prob in enumerate(probs):
                t_idx = valid_target_indices[i].item()
                if t_idx != q_idx:
                    results.append((prob.item(), t_idx, r_idx))

        results.sort(key=lambda x: x[0], reverse=True)

        for prob, t_idx, r_idx in results[:args.top_k]:
            t_name = g.vs[t_idx]["name"]
            rel_name = [r_idx]
            print(f"[{rel_name}] -> {t_name} | Score: {prob:.4f}")