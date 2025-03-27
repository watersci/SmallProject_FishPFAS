🔬 PFAS Tissue Cluster Analysis in New Jersey (Morone americana)

This small project explores PFAS (per- and polyfluoroalkyl substances) in fish tissue collected across New Jersey using open data from the EPA’s Water Quality Portal. Focusing on American white perch (Morone americana), the goal was to examine whether PFAS profiles have changed over time — and if so, how.

Using R and open tools, the project:

Downloads and processes PFAS data for fish tissue in New Jersey

Normalizes values to a common unit (pg/g) and handles non-detects using ½ the detection limit

Aggregates data by site, date, and compound and reshapes it into a matrix

Performs hierarchical clustering to identify groups of samples with similar PFAS fingerprints

Uses PERMANOVA (adonis) to test whether clusters reflect statistically distinct PFAS compositions

Explores whether cluster identity correlates with sampling date

Visualizes compositional shifts over time using PCA and stacked bar plots of relative concentrations

🧪 Key Finding:
Around 2011, PFAS profiles shifted significantly — not just in total concentration, but in the types of compounds detected. This suggests real-world chemical replacement or regulatory changes in PFAS use.

🛠️ Tools Used:
dataRetrieval for downloading WQP data

tidyverse for data wrangling

factoextra + ggfortify for clustering and PCA

vegan::adonis2 for multivariate testing

ggplot2 for temporal and compositional visualizations
