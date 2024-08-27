**Gene Expression and Patient Survival Analysis**

This project looks at how gene expression data can relate to patient survival times.
It focuses on identifying significant genes that might influence survival outcomes 
for different cancer types, such as Pancreatic Adenocarcinoma (PAAD). The project
involves analyzing gene expression data, clustering patients based on these data, 
and examining differences in survival times between these clusters.


Reading Gene Data: I take raw gene expression data and normalize it to make comparisons 
between patients easier.
Clustering: Then I group patients into clusters based on their gene expression profiles
to see if there are different survival outcomes.
Statistical Analysis: Then I use tests like ANOVA to see if these clusters are meaningfully
different in terms of survival, and we use linear regression to find out which genes are
important.


**Results**

After running the analysis, I identified key findings:

Patient and Gene Data:

Total Patients: 183
Patients with Recorded Death Data: 95
Total Genes Analyzed: 20,531
Genes after Filtering: 20,072
Clustering Analysis:

Number of Clusters: 3
Mean Survival Time by Cluster:
Cluster 0: 0.75 (standard deviation: 1.91, median: 0.28)
Cluster 1: 0.47 (standard deviation: 0.85, median: 0.37)
Cluster 2: -0.23 (standard deviation: 0.83, median: -0.38)
ANOVA Test:

P-Value: 0.003 (suggests statistically significant differences in survival times
between clusters)

Top Genes Identified:

Top Genes: ['C17orf105|284067', 'PPP1R3A|5506', 'MSH4|4438', 'COL9A1|1297', 'S100A5|6276',
    'TSPY3|728137', 'ISL2|64843', 'CBWD6|644019', 'FGF18|8817', 'FNDC8|54752']

Gene Weights: [0.73, 0.58, 0.57, 0.18, 0.17, 0.17, 0.17, 0.14, 0.12, 0.10]

These results indicate that certain genes have a strong association with patient 
survival and may serve as potential biomarkers for prognosis in cancer treatment.






