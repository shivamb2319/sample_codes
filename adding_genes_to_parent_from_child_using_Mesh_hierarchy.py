import networkx as nx

def add_genes_to_parents(concept_ids, genes, concept_to_genes, ontology):
    parents_to_update = list(ontology.predecessors(concept_ids))
    
    for parent in parents_to_update:
        if parent in concept_to_genes:
            concept_to_genes[parent] = list(set(concept_to_genes[parent] + genes))
        else:
            concept_to_genes[parent] = genes

    for parent in parents_to_update:
        add_genes_to_parents(parent, genes, concept_to_genes, ontology)

def get_concept_ids(concept_name, concept_to_concept_ids):
    return concept_to_concept_ids.get(concept_name, [])

def main():
    concept_id_to_genes = {}
    concept_to_concept_ids = {}
    concept_id_to_concept = {}
    ontology = nx.DiGraph()
    root_concept_ids = {}

    with open("Mesh_raw_file.txt") as f:
        for line in f:
            concept, concept_id = line.strip("\n").split(";")
            concept = concept.lower()
            concept_id_to_concept[concept_id] = concept
            concept_to_concept_ids.setdefault(concept, []).append(concept_id)
            inner_words = concept_id.split(".")
            if len(inner_words) == 1:
                root_concept_ids[concept] = inner_words[0]
            else:
                parent_id = ".".join(inner_words[:-1])
                ontology.add_edge(parent_id, concept_id)

    with open("disease_genes_file.txt") as f:
        for line in f:
            parts = line.strip("\n").split(':')
            concept_name = parts[0].strip().lower()
            gene_ids = [int(gene_id) for gene_id in parts[1].split(',')]

            concept_ids = get_concept_ids(concept_name, concept_to_concept_ids)

            if not concept_ids:
                print(f"Concept '{concept_name}' not found in the ontology graph. Skipping.")
                continue

            for concept_id in concept_ids:
                if concept_name in concept_id_to_genes:
                    concept_id_to_genes[concept_id] = list(set(concept_id_to_genes[concept_id] + gene_ids))
                else:
                    concept_id_to_genes[concept_id] = gene_ids


    concept_ids_copy = list(concept_id_to_genes.keys())

    for concept_id in concept_ids_copy:
        add_genes_to_parents(concept_id, concept_id_to_genes[concept_id], concept_id_to_genes, ontology)

    with open("final_genes.txt", "w") as output_file:
        for concept_id, genes in concept_id_to_genes.items():
            concept_name = concept_id_to_concept.get(concept_id, concept_id)
            unique_genes = sorted(set(genes))
            output_file.write(f"{concept_name}: {', '.join(map(str, unique_genes))}\n")

if __name__ == "__main__":
    main()


with open('final_genes.txt', 'r') as file:
    lines = file.readlines()

disease_genes = {}

for line in lines:
    parts = line.strip().split(':')
    if len(parts) == 2:
        disease_name = parts[0].strip()
        genes = parts[1].strip().split(', ')
        if disease_name in disease_genes:
            disease_genes[disease_name].extend(genes)
        else:
            disease_genes[disease_name] = genes

for disease_name in disease_genes:
    disease_genes[disease_name] = list(set(disease_genes[disease_name]))

with open('Output.txt', 'w') as file:
    for disease_name, genes in disease_genes.items():
        file.write(f'{disease_name}: {", ".join(genes)}\n')
