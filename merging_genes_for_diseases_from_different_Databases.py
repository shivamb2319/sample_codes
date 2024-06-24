# Read data from OMIM and store it in a dictionary

file1_data = {}

with open('Database1.txt', 'r') as file1:
    for line in file1:
        parts = line.strip().split(': ')
        if len(parts) == 2:
            disease_name = parts[0]
            genes = parts[1].split(', ')
            file1_data[disease_name] = set(genes)

# Read data from DisGeNEt and update the dictionary with unique genes

with open('Database2.txt', 'r') as file2:
    for line in file2:
        parts = line.strip().split(': ')
        if len(parts) == 2:
            disease_name = parts[0]
            genes = parts[1].split(', ')
            if disease_name in file1_data:
                file1_data[disease_name].update(genes)
            else:
                file1_data[disease_name] = set(genes)

# Create a new text file for the output

with open('Output.txt', 'w') as output_file:
    for disease_name, genes in file1_data.items():
        genes_str = ', '.join(sorted(genes))  # Sort genes and join them
        output_file.write(f"{disease_name}: {genes_str}\n")
      
