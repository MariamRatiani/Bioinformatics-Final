from sklearn.cluster import KMeans
from scipy.spatial.distance import euclidean
import numpy as np


class ClusterTightnessAnalyzer:
    """
    This class calculates the "tightness" of clusters for different numbers of clusters,
    providing insight into how well the gene expression data clusters together.
    """

    def __init__(self, gene_expression_data):
        """
        Initialize the analyzer with gene expression data.

        :param gene_expression_data: A 2D list or array where rows represent patients and columns represent genes.
        """
        self.gene_expression_data = gene_expression_data
        self.cluster_diameter_map = {}

    def analyze_cluster_diameters(self, cluster_sizes=None):
        """
        Calculate the average cluster diameters for a range of cluster sizes.

        :param cluster_sizes: A list of integers representing the number of clusters to try.
                              If None, the default range [2, 3, 4, 5, 6, 7, 8, 9, 10] is used.
        :return: A dictionary mapping the number of clusters to the average cluster diameter.

        	•	Key: The number of clusters used in the KMeans algorithm (e.g., 2, 3, 4, etc.).
            •	Value: The overall average diameter across all clusters when that specific
            number of clusters is used.
        """
        if cluster_sizes is None:
            cluster_sizes = [2, 3, 4, 5, 6, 7, 8, 9, 10]

        for num_clusters in cluster_sizes:
            # Perform KMeans clustering
            kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(self.gene_expression_data)

            cluster_distances = {}
            num_patients_by_cluster = {}
            for i in range(num_clusters):
                cluster_distances[i] = 0
                num_patients_by_cluster[i] = 0

            # Calculate the sum of distances for each cluster
            for i in range(len(kmeans.labels_)):
                patient_cluster_label = kmeans.labels_[i]
                num_patients_by_cluster[patient_cluster_label] += 1
                patient_cluster_center = kmeans.cluster_centers_[patient_cluster_label]
                cluster_distances[patient_cluster_label] += euclidean(self.gene_expression_data[i],
                                                                      patient_cluster_center)

            # Compute the average distance (diameter) for each cluster
            for i in range(num_clusters):
                cluster_distances[i] /= float(num_patients_by_cluster[i])
            self.cluster_diameter_map[num_clusters] = sum(cluster_distances.values()) / float(num_clusters)

        return self.cluster_diameter_map

    def get_cluster_diameters(self):
        """
        Retrieve the calculated cluster diameters.

        :return: A dictionary mapping the number of clusters to the average cluster diameter.
        """
        return self.cluster_diameter_map
