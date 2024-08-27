import numpy as np
from sklearn.cluster import KMeans
from sklearn import linear_model
from scipy import stats

class ClusterAnalysis:
    def __init__(self, num_clusters, patient_ids, patient_death_data, gene_expression_data, gene_labels):
        self.num_clusters = num_clusters
        self.patient_ids = patient_ids
        self.patient_death_data = patient_death_data
        self.gene_expression_data = gene_expression_data
        self.gene_labels = gene_labels
        self.kmeans_model = None
        self.deaths_by_cluster = {}
        self.gene_expression_by_cluster = {}
        self.cluster_means = {}
        self.cluster_stddevs = {}
        self.cluster_medians = {}
        self.anova_p_value = None
        self.top_genes = []
        self.top_gene_weights = []

    def run_kmeans(self):
        """Runs KMeans clustering and initializes the cluster data structures."""
        print("RUN KMEANS WITH N CLUSTERS: ", self.num_clusters)
        self.kmeans_model = KMeans(n_clusters=self.num_clusters, random_state=0).fit(self.gene_expression_data)
        for i in range(self.num_clusters):
            self.deaths_by_cluster[i] = []
            self.gene_expression_by_cluster[i] = []
        self.assign_patients_to_clusters()

    def assign_patients_to_clusters(self):
        """Assigns patients and their gene expression data to the respective clusters."""
        print("SPLIT PATIENTS ON CLUSTER")
        for i in range(len(self.kmeans_model.labels_)):
            cluster_label = self.kmeans_model.labels_[i]
            patient_id = self.patient_ids[i]
            if patient_id in self.patient_death_data.keys():
                self.deaths_by_cluster[cluster_label].append(self.patient_death_data[patient_id])
                self.gene_expression_by_cluster[cluster_label].append(self.gene_expression_data[i])

    def calculate_statistics(self):
        """Calculates mean, median, and standard deviation for each cluster."""
        for i in range(self.num_clusters):
            self.cluster_means[i] = np.mean(np.array(self.deaths_by_cluster[i]))
            self.cluster_stddevs[i] = np.std(np.array(self.deaths_by_cluster[i]))
            self.cluster_medians[i] = np.median(np.array(self.deaths_by_cluster[i]))

        print("mean: ", self.cluster_means)
        print("standard deviation: ", self.cluster_stddevs)
        print("median: ", self.cluster_medians)

    def perform_anova(self):
        """Performs ANOVA to test the statistical significance of differences between clusters."""
        print("PERFORMING ANOVA")
        if self.num_clusters > 2:
            self.anova_p_value = stats.f_oneway(*[self.deaths_by_cluster[i] for i in range(self.num_clusters)]).pvalue
        print("ANOVA P VALUE: ", self.anova_p_value)
        return self.anova_p_value

    def analyze_clusters(self):
        """Main method to run analysis if ANOVA shows significant differences between clusters."""
        if self.anova_p_value < 0.05:
            max_mean_cluster = max(self.cluster_means, key=self.cluster_means.get)

            label_data, gene_data = self.create_label_and_gene_data(max_mean_cluster)
            self.run_linear_regression(gene_data, label_data)
            self.compare_gene_expression_in_clusters(max_mean_cluster)

    def create_label_and_gene_data(self, max_mean_cluster):
        """Creates label data for regression where patients in max_mean_cluster are labeled 1, others 0."""
        label_data = []
        gene_data = []
        for i in range(self.num_clusters):
            for j in range(len(self.gene_expression_by_cluster[i])):
                label_data.append(1 if i == max_mean_cluster else 0)
                gene_data.append(self.gene_expression_by_cluster[i][j])
        return label_data, gene_data

    def run_linear_regression(self, gene_data, label_data):
        """Performs linear regression to find the most important genes for cluster membership."""
        reg_model = linear_model.LinearRegression()
        reg_model.fit(gene_data, label_data)
        absolute_coefficients = list(map(abs, reg_model.coef_.tolist()))

        # Find the top 10 most influential genes
        for i in range(10):
            max_ind = absolute_coefficients.index(max(absolute_coefficients))
            self.top_genes.append(self.gene_labels[max_ind])
            self.top_gene_weights.append(max(absolute_coefficients))
            absolute_coefficients[max_ind] = 0

    def compare_gene_expression_in_clusters(self, max_mean_cluster):
        """Compares gene expression of the most influential gene between the max cluster and others."""
        top_gene_index = self.gene_labels.index(self.top_genes[0])
        max_cluster_expression = 0
        other_clusters_expression = 0
        # max_cluster_size = len(self.gene_expression_by_cluster[max_mean_cluster])

        for i in range(self.num_clusters):
            for j in range(len(self.gene_expression_by_cluster[i])):
                if i == max_mean_cluster:
                    max_cluster_expression += self.gene_expression_by_cluster[i][j][top_gene_index]
                else:
                    other_clusters_expression += self.gene_expression_by_cluster[i][j][top_gene_index]

    def get_top_genes(self):
        """Returns the top 10 genes and their weights identified by the analysis."""
        return self.top_genes, self.top_gene_weights