# cancer_prostate

## Multiomic Data Analytics Integration in Prostate Cancer

### Background  

Prostate cancer (PCa), a hormone-dependent oncological disease, is a deadly disease for men associated with heterogeneous clinical outcome. Clinically, it is characterized by indolent phenotypes, rapid progression and aggressive metastatic disease (Sathianathen et al., 2018). Among cancer worldwide, PCa as a common male tumor has an incidence of second lung cancer (Siegel et al., 2019). This neoplastic process according to the Global Cancer Incidence, Mortality and Prevalence (GLOBOCAN) has become a leading cause of cancer worldwide, with over 1.4 million new cases and a total of 375,000 deaths in 2020 (Sung et al., 2020). The heterogeneity of the disease calls for the use of all types of omics data necessary to promote precision medicine. PCa progression is slow compared to other tumors but still harms patient long-term health (Eggener et al., 2015). The multifocal and heterogeneity of primary PCa is associated with poor prognosis (wang et al., 2018), this means that intervention for patients with metastasis are still urgently needed. While patients with slow disease progression show better outcomes, patients with higher aggressiveness show worse treatment prognosis and current clinicopathological indicators do not distinguish well between patients based on their outcomes at the initial stage of the disease (Zhang et al., 2020). Consequently, efficient identification of the risk level of prostate cancer patients and precise therapeutic targets has always been an important equation to solve in PCa. The occurrence and development of PCa is related to multisystem and multilevel pathological changes. Studies at single omics level often have limitations, while combined analysis of multiple omics data can better and more comprehensively develop targeted markers for PCa therapy. A study done by Sinha et al PCa multi omics revealed that proteomics features were significantly more informative for biochemical recurrence prediction than genomic, epigenomics or transcriptomics (Sinha et al., 2019). As much as the single technologies have identified and shed light on the mechanisms of the tumor progression, subtypes and finding new treatment targets, a holistic molecular perspective of the mechanisms, prospective biomarkers and drug targets and responses can only be uncovered when a systems biology approach is adopted. Hence in the current study we shall be using good quality multi-omics data sets from public databases to gain a better understanding of tumor progression, subtyping and finding novel biomarkers that potentially address individual variations in drug responses among prostate cancer patients. In addition, we shall also provide a simplified protocol for routine integration of multi-omics data sets to answer biological questions.


### Aim
Identify African descent biomarkers associated with Prostate cancer using multi omics available data for better management of patients. 

### Objectives

1.Explain prostate cancer biological complexity using multi omics approach

2.Provide a simplified and standardized African specific prostate cancer management 

### Outcome Measures 

1.	Genetic profiling results - Mutation identified via whole-genome sequencing will be recorded.
2.	Transcriptional profiling results - Determining the transcriptomic information of prostate cancer
3.	Epigenomic profiling results - Determining the epigenomic status of specific genes
4.	Proteomic profiling results – determine the proteomic information of the prostate cancer
5.	Metabolomic profiling results – determine the profiles of the metabolites in prostate cancer patients

### Methods

#### Workflow

-Single omics-

![Single omics pathways](https://user-images.githubusercontent.com/102041566/161955431-539e2ed4-c488-4e74-978e-ab74649aa6aa.PNG)


- Multi omics -

![image](https://user-images.githubusercontent.com/102041566/161955384-bdb603ec-f4fc-4853-9ddb-937009b3a510.png)


### Tools and software

Whole genome sequence: FastQC, BWA, Bcftools, snpEff 

RNASeq:feature count, HTSeq, Sleuth, 

Metabolome: [MetaboAnalytR](https://www.metaboanalyst.ca/docs/RTutorial.xhtml)

Single-omics analysis will be done on Linux and multi-omics analysis will be done in R. 

### Input and output

Output from single omics will be .tsv format based on different features(SNPs, transcripts, metabolome)



## Team structure

### Core members

1-Team lead: Zedias Chikwambi

2-Co-lead: [Marie Hidjo](https://github.com/mariehidjo/cancer_prostate)

3-Technical lead: [Pageneck Chikondowa](github.com/pageneck)

4-Writer (Github): [Glory Jayeoba](https://github.com/gloryife)

5-Writer (Manuscript): [Lawrence Afolabi](https://github.com/itslawrenceb)

### Scripting and Workflow team 

1-[David enoma](https://github.com/davidenoma)

2-Glory Jayeoba

3-[Marie Hidjo](https://github.com/mariehidjo/cancer_prostate)

4-[Vincent Aketch Nyangwara](https://github.com/vinaketch)

5-Pageneck Chikondowa

### Manuscript team 

1-Zedias Chikwambi

2-Lawrence Afolabi

3-David Juma

### Slides 

1-Zedias Chikwambi

2-Marie Hidjo



### References 

C. Manzoni, D. A. Kia, J. Vandrovcova et al., “Genome, transcriptome and proteome: the rise of omics data and their integration in biomedical sciences,” Briefings in Bioinformatics, vol. 19, no. 2, pp. 286–302, 2018. 

Sung, H.; Ferlay, J.; Siegel, R.L.; Laversanne, M.; Soerjomataram, I.; Jemal, A.; Bray, F. Global Cancer Statistics 2020: GLOBOCAN Estimates of Incidence and Mortality Worldwide for 36 Cancers in 185 Countries. CA Cancer J. Clin. 2021, 71, 209–249. 

Sathianathen, N.J.; Konety, B.R.; Crook, J.; Saad, F.; Lawrentschuk, N. Landmarks in Prostate Cancer. Nat. Rev. Urol. 2018, 15, 627–642. 

Gómez-Cebrián, N.; Poveda, J.L.; Pineda-Lucena, A.; Puchades-Carrasco, L. Metabolic Phenotyping in Prostate Cancer Using Multi-Omics Approaches. Cancers 2022, 14, 596. https://doi.org/10.3390/ cancers14030596
R.L. Siegel, K.D. Miller, A. Jemal, Cancer statistics, 2019, CA Cancer J. Clin. 69 (1)
(2019), https://doi.org/10.3322/caac.21551 (PubMed PMID: 30620402). 

P. Mertins, D.R. Mani, K.V. Ruggles, M.A. Gillette, K.R. Clauser, P. Wang. Proteogenomics connects somatic mutations to signalling in breast cancer, Nature 534 (7605) (2016) 55–62, https://doi.org/10.1038/nature18003 (PubMed PMID: 27251275).

 H. Zhang, T. Liu, Z. Zhang, S.H. Payne, B. Zhang, J.E. McDermott, et al., Integrated
proteogenomic characterization of human high-grade serous ovarian cancer, Cell 166 (3) (2016) 755–765, https://doi.org/10.1016/j.cell.2016.05.069 (PubMed PMID: 27372738).

S.E. Eggener, A.S. Cifu, C. Nabhan, Prostate cancer screening, JAMA 314 (8) (2015) 825–826, https://doi.org/10.1001/jama.2015.8033 (PubMed PMID:26305653). 

G. Wang, D. Zhao, D.J. Spring, R.A. DePinho, Genetics and biology of prostate cancer, Genes Dev. 32 (17–18) (2018) 1105–1140, https://doi.org/10.1101/gad. 315739.118 (PubMed PMID: 30181359). 
R.L. Siegel, K.D. Miller, A. Jemal, Cancer statistics, 2019, CA Cancer J. Clin. 69 (1) (2019), https://doi.org/10.3322/caac.21551 (PubMed PMID: 30620402). 

A. Sinha, V. Huang, J. Livingstone, J. Wang, N.S. Fox, N. Kurganovs, et al., The Proteogenomic landscape of curable prostate cancer, Cancer Cell 35 (3) (2019), https://doi.org/10.1016/j.ccell.2019.02.005 (PubMed PMID: 30889379). 



