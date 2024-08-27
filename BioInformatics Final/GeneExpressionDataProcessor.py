import csv
import numpy as np


# rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data
# row - gene_id
# column - patient


class GeneExpressionDataProcessor:
    """
    This class handles the processing of gene expression data for patients with a specific type of cancer.
    It provides functionality to read, filter, and normalize gene expression data, as well as to extract
    gene expression data for selected genes.
    """

    def __init__(self):
        self.all_gene_expression = []

    def get_patients_genes_and_gene_expression(self, cancer, patient_death_data):
        """
        Reads the gene expression data for patients with a specified cancer type.
        Filters the data to include only patients who have recorded days-to-death information.
        Normalizes the gene expression values for comparison purposes.

        :param cancer: The type of cancer (e.g., "PAAD" for Pancreatic Adenocarcinoma).
        :param patient_death_data: A dictionary containing patient IDs and their days-to-death values.
        :return: A tuple containing filtered patient IDs, gene names, and the normalized gene expression data.
        """
        line_num = 0
        gene_expression_all = []
        genes = []
        patients = []

        # Lists to store the maximum and minimum values for each gene across all patients
        # These will be used for normalization purposes
        max_val_genes = []
        min_val_genes = []

        print("READING THE DATA")
        # Open the gene expression file for the specified cancer type
        with open(
                cancer + ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt") as tsv:
            for line in csv.reader(tsv, dialect="excel-tab"):
                if line_num == 0:
                    # The first row contains patient IDs
                    patients = line[1:]
                    line_num += 1
                elif line_num == 1:
                    # Skip the second row (usually contains headers or meta-data)
                    line_num += 1
                else:
                    # Convert gene expression values from strings to floats
                    gene_data = [float(d) for d in line[1:]]
                    # Exclude genes that have no expression across all patients
                    if max(gene_data) != 0:
                        genes.append(line[0])
                        gene_expression_all.append(gene_data)
                        # Store max and min values for normalization
                        max_val_genes.append(max(gene_data))
                        min_val_genes.append(min(gene_data))
                    line_num += 1

        # Extract the unique patient label from the patient ID (usually part of the ID string)
        patients = [patient.split("-")[2].lower() for patient in patients]
        patients_filtered = []

        print("FILTER ONLY PATIENTS WITH DEATH")
        # Transpose gene expression data to have patients as rows and genes as columns
        gene_expression_transposed = np.asarray(gene_expression_all).T.tolist()
        gene_expression_transposed_filtered = []

        # Filter to include only patients with recorded death data
        for i in range(len(patients)):
            if patients[i] in patient_death_data:
                patients_filtered.append(patients[i])
                gene_expression_transposed_filtered.append(gene_expression_transposed[i])

        # Transpose back to have genes as rows and patients as columns
        gene_expression = np.asarray(gene_expression_transposed_filtered).T.tolist()

        print("NORMALIZE THE DATA")
        # Normalize gene expression values to a range of 0 to 1 for each gene
        for a in range(len(gene_expression)):
            for c in range(len(gene_expression[a])):
                gene_expression[a][c] = (gene_expression[a][c] - min_val_genes[a]) / (
                        max_val_genes[a] - min_val_genes[a])

        print("NUM PATIENTS ", len(patients))
        print("NUM PATIENTS FILTERED", len(patients_filtered))
        print("NUM GENES BEFORE " + str(line_num - 2))
        print("NUM GENES FILTERED", len(gene_expression))

        self.all_gene_expression = gene_expression

        return patients_filtered, genes, gene_expression

    def get_gene_expression_for_highest_variation_genes(self, all_gene_expression_data, selected_gene_names,
                                                        all_gene_names):
        """
        Extracts gene expression data for the genes with the highest variation, as specified in selected_gene_names.

        :param all_gene_expression_data: The complete gene expression dataset for all genes and patients.
        :param selected_gene_names: A list of gene names that have the highest variation.
        :param all_gene_names: A list of all gene names present in the dataset.
        :return: A subset of the gene expression data containing only the selected genes.
        """
        selected_gene_indices = []

        # Collect indices of selected genes from the full gene list
        for gene_name in selected_gene_names:
            selected_gene_indices.append(all_gene_names.index(gene_name))

        print("Selected Gene Indices: ", selected_gene_indices)

        # Transpose the gene expression data to have patients as rows and genes as columns
        transposed_gene_expression_data = np.asarray(all_gene_expression_data).T.tolist()

        # Filter and collect the gene expression data for the selected genes only
        selected_gene_expression_data = []
        for patient_data in transposed_gene_expression_data:
            selected_gene_expression_data.append([patient_data[index] for index in selected_gene_indices])

        return selected_gene_expression_data