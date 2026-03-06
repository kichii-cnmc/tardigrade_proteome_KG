import torch
import torch.nn.functional as F
import igraph as ig
import numpy as np
import argparse
from torch_geometric.nn import RGCNConv

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