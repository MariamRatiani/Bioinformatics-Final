import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.mlab import PCA
import numpy as np
from sklearn import datasets, linear_model

from ClusterAnalysis import ClusterAnalysis
from ClusterDiameterPlotter import ClusterDiameterPlotter
from ClusterTightnessAnalyzer import ClusterTightnessAnalyzer
from GeneExpressionDataProcessor import GeneExpressionDataProcessor
from PatientDeathDataProcessor import PatientDeathDataProcessor

"""
    ორი დათა ფაილი გვაქვს: 1 -> პაციენტები რამდენ დღეში მოკვდნენ
                           2 -> პაციენტები და მათში გენების ექსპრესიის რიცხვები 
"""

cancer = "PAAD"


# goal is to highlight the top genes that are most strongly associated
# with patient survival times
def regression_death_data(patients_to_death, patients, genes, gene_expr):
    gene_transp = np.asarray(gene_expr).T.tolist()

    # standardized days-to-death values for each patient
    label_data = []
    # the gene expression data for each patient
    gene_data = []
    for i in range(len(patients)):
        patient_id = patients[i]
        label_data = label_data + [patients_to_death[patient_id]]
        gene_data.append(gene_transp[i])

    reg_model = linear_model.LinearRegression()
    reg_model.fit(gene_data, label_data)
    model_coef_abs_map = list(map(abs, reg_model.coef_.tolist()))
    max_weights_gene_ind = []
    max_weights_gene = []
    for i in range(20):
        max_ind = model_coef_abs_map.index(max(model_coef_abs_map))
        max_weights_gene_ind.append(max_ind)
        max_weights_gene.append(genes[max_ind])
        model_coef_abs_map[max_ind] = 0
    print("my output important genes:", max_weights_gene)
    return max_weights_gene


def main():
    patient_death_data_processor = PatientDeathDataProcessor(cancer)
    gene_data_processor = GeneExpressionDataProcessor()

    # dictionary of patient ids and their death data
    patient_death_data = patient_death_data_processor.get_standardized_patients_to_death()
    patients, genes, gene_expression = gene_data_processor.get_patients_genes_and_gene_expression(cancer,
                                                                                                  patient_death_data)

    max_weights_gene = regression_death_data(patient_death_data, patients, genes, gene_expression)

    # patients -> პაციენტები, რომლებზეც გვაქვს ინფო სიკვდილისა და გენების ექსპრესიის შესახებ

    gene_names = max_weights_gene
    gene_expr_max_var = gene_data_processor.get_gene_expression_for_highest_variation_genes(gene_expression, gene_names,
                                                                                            genes)

    # num_clusters_to_diameter = cluster_diameters_by_size(gene_expr_max_var)
    cluster_analyzer = ClusterTightnessAnalyzer(gene_expr_max_var)
    num_clusters_to_diameter = cluster_analyzer.analyze_cluster_diameters()

    plotter = ClusterDiameterPlotter(num_clusters_to_diameter, cancer + " Death")
    plotter.plot()
    # plot_diameter_vs_num_clusters(num_clusters_to_diameter, cancer + " Death")

    num_clusters = 3
    # get_cluster_stats(num_clusters, patients, patient_death_data, gene_expr_max_var, gene_names, genes)
    cluster_analysis = ClusterAnalysis(num_clusters, patients, patient_death_data, gene_expr_max_var, gene_names)
    # Run the KMeans clustering
    cluster_analysis.run_kmeans()
    # Calculate statistical measures for the clusters
    cluster_analysis.calculate_statistics()
    # Perform ANOVA to check for significant differences between clusters
    if cluster_analysis.perform_anova() < 0.05:
        # If significant, analyze clusters to find top genes
        cluster_analysis.analyze_clusters()

        # Get the top genes that are most associated with cluster membership
        top_genes, top_gene_weights = cluster_analysis.get_top_genes()
        print("Top Genes: ", top_genes)
        print("Top Gene Weights: ", top_gene_weights)


main()

