import matplotlib.pyplot as plt

# DONE

class ClusterDiameterPlotter:
    def __init__(self, num_clusters_to_diameter, cancer_type):
        """
        Initializes the ClusterDiameterPlotter with the necessary data.

        :param num_clusters_to_diameter: Dictionary where keys are the number of clusters and values are the corresponding average diameters.
        :param cancer_type: The type of cancer being analyzed (used for plot title).
        """
        self.num_clusters_to_diameter = num_clusters_to_diameter
        self.cancer_type = cancer_type

    def __putLabels(self):
        plt.ylabel("Diameter")
        plt.xlabel("Number of Clusters")
        plt.title("Average Diameter Across Cluster Sizes " + self.cancer_type)


    def plot(self):
        """
        Plots the average cluster diameter against the number of clusters.
        """
        x_is = list(self.num_clusters_to_diameter.keys())
        y_is = [self.num_clusters_to_diameter[cluster_size] for cluster_size in x_is]

        plt.plot(x_is, y_is)
        self.__putLabels()

        plt.axis([0, 10, 0, max(self.num_clusters_to_diameter.values())])
        plt.grid(True)
        plt.show()
