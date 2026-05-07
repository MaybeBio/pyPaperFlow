- [PMID 36851914 - Computational prediction of disordered binding regions.](#pmid-36851914---computational-prediction-of-disordered-binding-regions)
  - [Title](#title)
  - [Pmid](#pmid)
  - [Keywords](#keywords)
  - [Abstract](#abstract)
  - [Conclusion](#conclusion)
- [PMID 38701796 - DR-BERT: A protein language model to annotate disordered regions.](#pmid-38701796---dr-bert-a-protein-language-model-to-annotate-disordered-regions)
  - [Title](#title-1)
  - [Pmid](#pmid-1)
  - [Keywords](#keywords-1)
  - [Abstract](#abstract-1)
- [PMID 39763873 - A deep learning method for predicting interactions for intrinsically disordered regions of proteins.](#pmid-39763873---a-deep-learning-method-for-predicting-interactions-for-intrinsically-disordered-regions-of-proteins)
  - [Title](#title-2)
  - [Pmid](#pmid-2)
  - [Keywords](#keywords-2)
  - [Abstract](#abstract-2)
  - [Discussion](#discussion)
    - [Uses and limitations](#uses-and-limitations)
    - [Challenges](#challenges)
    - [Performance of AF2 and AF3](#performance-of-af2-and-af3)
    - [Future Directions](#future-directions)
- [PMID 40286477 - Amino acid sequence-based IDR classification using ensemble machine learning and quantum neural networks.](#pmid-40286477---amino-acid-sequence-based-idr-classification-using-ensemble-machine-learning-and-quantum-neural-networks)
  - [Title](#title-3)
  - [Pmid](#pmid-3)
  - [Keywords](#keywords-3)
  - [Abstract](#abstract-3)
- [PMID 41378882 - Impact of intrinsically disordered regions and functional disorder hotspots in the human kinome.](#pmid-41378882---impact-of-intrinsically-disordered-regions-and-functional-disorder-hotspots-in-the-human-kinome)
  - [Title](#title-4)
  - [Pmid](#pmid-4)
  - [Keywords](#keywords-4)
  - [Abstract](#abstract-4)
  - [Conclusion](#conclusion-1)
- [PMID 41534519 - Disobind: A sequence-based, partner-dependent contact map and interface residue predictor for intrinsically disordered regions.](#pmid-41534519---disobind-a-sequence-based-partner-dependent-contact-map-and-interface-residue-predictor-for-intrinsically-disordered-regions)
  - [Title](#title-5)
  - [Pmid](#pmid-5)
  - [Keywords](#keywords-5)
  - [Abstract](#abstract-5)

---

<a id="pmid-36851914-computational-prediction-of-disordered-binding-regions"></a>

# PMID 36851914 - Computational prediction of disordered binding regions.

## Title

Computational prediction of disordered binding regions.

## Pmid

36851914

## Keywords

CAID, Critical Assessment of Intrinsic Disorder, CASP, Critical Assessment of techniques for protein Structure Prediction, DL, deep learning, Disordered binding regions, IDP, intrinsically disordered protein, IDR, intrinsically disordered region, Intrinsic disorder, ML, machine learning, MoRF, molecular recognition fragment, Molecular recognition features, NN, neural network, Protein-lipid interactions, Protein-nucleic acids interactions, Protein-protein interactions, SLiM, short linear sequence motif, Short linear motifs

## Abstract

One of the key features of intrinsically disordered regions (IDRs) is their ability to interact with a broad range of partner molecules. Multiple types of interacting IDRs were identified including molecular recognition fragments (MoRFs), short linear sequence motifs (SLiMs), and protein-, nucleic acids- and lipid-binding regions. Prediction of binding IDRs in protein sequences is gaining momentum in recent years. We survey 38 predictors of binding IDRs that target interactions with a diverse set of partners, such as peptides, proteins, RNA, DNA and lipids. We offer a historical perspective and highlight key events that fueled efforts to develop these methods. These tools rely on a diverse range of predictive architectures that include scoring functions, regular expressions, traditional and deep machine learning and meta-models. Recent efforts focus on the development of deep neural network-based architectures and extending coverage to RNA, DNA and lipid-binding IDRs. We analyze availability of these methods and show that providing implementations and webservers results in much higher rates of citations/use. We also make several recommendations to take advantage of modern deep network architectures, develop tools that bundle predictions of multiple and different types of binding IDRs, and work on algorithms that model structures of the resulting complexes.

## Conclusion

IDRs interact with many different molecular partners including proteins, DNA, RNA, lipids, small molecules, carbohydrates, and metals. The knowledge of these interactions is rather limited, which motivates development of computations tools that predicts them from the readily available protein sequences. This comprehensive survey of sequence-based predictors of binding IDRs covers a wide range of interacting partners. We identify and summarize a large collection of 38 predictors that consider 5 different types of interacting IDRs. The MoRF predictors are the largest category with 21 methods, followed by 9 SLiM predictors, 3 predictors of protein-binding IDRs, 3 methods that predict protein/DNA/RNA binding IDRs and 2 predictors of lipid-binding IDRs. We find that these methods rely on a diverse range of predictive architectures that include scoring functions, regular expressions, machine learning models and meta-predictors, where about three-quarters of them utilize machine learning algorithms. We observe a couple of recent trends to develop deep network-based models and to extend coverage to new types of interacting IDRs, such as RNA, DNA and lipid binding regions. We also note a high rate of availability of these methods, with over 70% that are provided to the end users as either webservers and/or standalone code. Furthermore, we analyze relation between availability and impact/use of these methods. We find that methods which are more broadly available, as both webserver and source code, are substantially more cited/used when compared to those that are available in either format, while methods that do not offer a publicly available implementation suffer low use/citations. Moreover, we also find that the availability should be maintained since tools that were originally made available and are currently not functional observe a large drop in the use/citations. The latter observations strongly suggest that future predictors should be made available in both formats upon publication and should be maintained after publication.

While IDPs interact with a broad range of molecular partners, we show that the current predictors are largely focused on two types of binding IDRs, MoRFs and SLiMs. A particularly acute situation concerns prediction of nucleic acid and lipid-binding IDRs, where only a handful of methods are available. The prediction of small molecules-, carbohydrates-, and metal-binding IDRs is not feasible at the moment, given a very small amount of ground truth data. The need to develop new predictors of DNA and RNA binding regions is further motivated by the inclusion of this prediction category in the pending CAID2 experiment. Consequently, one of the key future directions would be to diversify the development efforts to more uniformly cover different types of binding IDRs.

Results of the recently completed CAID assessment show that predictors of binding IDRs offer modest levels of predictive performance [63], suggesting that there is a large room for improvement. We observe that none of the methods that participated in this evaluation use deep learning models. The recent influx of the deep learning-based predictors of binding IDRs will likely result in improved predictive quality. This claim stems from a recent study that empirically demonstrates that deep learning-based predictors of intrinsic disorder significantly outperform other types of models [64]. The drive to use deep learning models is also motivated by the growing and successful use of these models in related areas of bioinformatics [148], such as prediction of protein-protein interactions 149, 150, 151 and protein function 152, 153, 154. We envision that majority of future predictors of binding IDRs will likely rely on deep neural networks. We encourage the developers to consider modern network topologies, such as the recently developed transformers [155], that were used to very accurately predict protein structures [156].

Some IDPs include IDRs that interact with different types of ligands and yet most of the current methods cover a single ligand type. Consequently, users are forced to use multiple methods and convert between different output formats to obtain a complete prediction. These difficulties could be alleviated with solutions that bundle multiple predictors, however, the only such solution to date is the DEPICTER webserver [157]. Moreover, there are only a handful of methods that predict IDRs that bind to multiple ligand types, such as DisoRDPbind, flDPnn and DeepDISObind, that target protein, RNA and DNA-binding IDRs. Consequently, we advocate for the development of new tools that address predictions of multiple and many different types of binding IDRs. Furthermore, some IDRs can bind multiple partner types, which corresponds to multi-label (multi-output) learning. Prediction of such multifunctional IDRs is possible with the DMRpred method, although this tool does not provide types of binding partners [158]. Thus, new tools that would cast this prediction as multi-labels problem should be developed. We note that multi-labels predictors are widely used in related areas, such as prediction of subcellular localization 159, 160, 161, 162, nucleic acid binding proteins [163], enzymatic functions [164], and ion channel types [165].

Prediction of the binding IDRs in protein sequences should be followed by modelling structures of the corresponding complexes (i.e., IDRs fold upon binding). While computational protein docking has been extensively pursued over the past several years [166], studies that investigate docking with IDPs are lagging behind since IDPs are difficult to model. Daisuke Kihara’s lab developed a pioneering approach for IDP-protein docking, IDP-LZerD 167, 168. This method produces a docking model from the 3D structure of the receptor and the sequence of interacting IDP. Docking an IDP is conceptually similar to protein-small peptide docking, but technically more challenging because conformation of the IDP on the receptor’s surface has to be predicted. In IDP-LZerD, this is done by docking and stitching short protein fragments taken from the binding IDR. Moreover, a recent benchmark study that evaluates three methods capable of docking with IDPs, IDP-LZerD 167, 168, CABS-Dock [169] and AlphaFold-Multimer [170], shows that they accurately identify location of the binding site but struggle with atomic-levels details of the structure [171], suggesting that further research is needed.

Lastly, databases like D2P2
[172], MobiDB 173, 174, 175, 176 and DescribePROT [177] provide convenient access to pre-computed predictions of disorder for millions of proteins. However, they typically contain a limited number of binding IDR predictions, with DescribePROT covering the most diverse range that includes putative protein, RNA and DNA-binding IDRs. This coverage should be extended in the future as more methods that cover a broader range of binding IDRs will be developed. In turn, this effort motivates the development of runtime-efficient predictors that can be used to perform predictions on such large scale. Examples of current fast tools include ANCHOR2, DisoRDPbind and fMoRFpred, that were shown to produce predictions in about 1 second per protein in the CAID experiment [63].


<!-- PAPER_BREAK -->

---

<a id="pmid-38701796-dr-bert-a-protein-language-model-to-annotate-disordered-regions"></a>

# PMID 38701796 - DR-BERT: A protein language model to annotate disordered regions.

## Title

DR-BERT: A protein language model to annotate disordered regions.

## Pmid

38701796

## Keywords

IDP, IDR, deep learning, disorder, machine learning, protein language model, protein structure prediction

## Abstract

Despite their lack of a rigid structure, intrinsically disordered regions (IDRs) in proteins play important roles in cellular functions, including mediating protein-protein interactions. Therefore, it is important to computationally annotate IDRs with high accuracy. In this study, we present Disordered Region prediction using Bidirectional Encoder Representations from Transformers (DR-BERT), a compact protein language model. Unlike most popular tools, DR-BERT is pretrained on unannotated proteins and trained to predict IDRs without relying on explicit evolutionary or biophysical data. Despite this, DR-BERT demonstrates significant improvement over existing methods on the Critical Assessment of protein Intrinsic Disorder (CAID) evaluation dataset and outperforms competitors on two out of four test cases in the CAID 2 dataset, while maintaining competitiveness in the others. This performance is due to the information learned during pretraining and DR-BERT's ability to use contextual information.


<!-- PAPER_BREAK -->

---

<a id="pmid-39763873-a-deep-learning-method-for-predicting-interactions-for-intrinsically-disordered-regions-of-proteins"></a>

# PMID 39763873 - A deep learning method for predicting interactions for intrinsically disordered regions of proteins.

## Title

A deep learning method for predicting interactions for intrinsically disordered regions of proteins.

## Pmid

39763873

## Keywords

Intrinsically disordered proteins (IDP), deep learning (DL), intrinsically disordered regions (IDR), protein language model (pLMs), protein structure

## Abstract

Intrinsically disordered proteins or regions (IDPs/IDRs) adopt diverse binding modes with different partners, from coupled-folding-and-binding, to fuzzy binding, to fully-disordered binding. Characterizing IDR interfaces is challenging experimentally and computationally. The state-of-the-art AlphaFold-multimer and AlphaFold3 can be used to predict IDR binding sites, although they are less accurate at their benchmarked confidence cutoffs. Here, we developed Disobind, a deep-learning method that predicts inter-protein contact maps and interface residues for an IDR and its partner, given their sequences. It uses sequence embeddings from the ProtT5 protein language model. Disobind outperforms state-of-the-art interface predictors for IDRs. It also outperforms AlphaFold-multimer and AlphaFold3 at multiple confidence cutoffs. Combining Disobind and AlphaFold-multimer predictions further improves the performance. In contrast to current methods, Disobind considers the context of the binding partner and does not depend on structures and multiple sequence alignments. Its predictions can be used to localize IDRs in large assemblies and characterize IDR-mediated interactions.

## Discussion

Here, we discuss about the uses and limitations of Disobind, challenges involved, the performance of AlphaFold, and future directions.


### Uses and limitations

  Predictions from Disobind+AF2 may be used to improve the localization of IDRs in integrative models of large assemblies, our primary motivation for developing the method. Macromolecular assemblies contain significant portions of IDRs, for example, the Fg Nups in the nuclear pore complex, the MBD3-IDR in the nucleosome remodeling and deacetylase (NuRD) complex, and the N-terminus of Plakophilin1 in the desmosome (Akey et al., 2022; Arvindekar et al., 2022; Pasani et al., 2024). These regions typically lack data for structural modeling, resulting in integrative models of poor precision (Arvindekar et al., 2022; Pasani et al., 2024). The predictions from Disobind+AF2 can be used in integrative modeling methods such as IMP, HADDOCK, and Assembline as inter-protein distance restraints (Dominguez et al., 2003; Rantos et al., 2022; Russel et al., 2012). Even the coarse-grained contact map and interface residue predictions would be useful in such cases to improve the precision of these regions in the integrative model. Further, the predicted contacts can be combined with molecular dynamics (MD) simulations to generate ensembles of IDRs in complexes, providing mechanistic insights into their dynamic behaviour. Additionally, our method can be used to characterize interactions involving IDRs across proteomes. This may aid in identifying new binding motifs for IDRs, potentially linked to their sub-cellular localization or function. Finally, predictions from our method may aid in modulating interactions involving IDRs, for example, by suggesting plausible mutations.

  However, our method also has several limitations. First, it is limited to binary IDR-partner complexes. For complexes formed by three or more subunits, binary predictions from Disobind would need to be combined. Second, it assumes that the IDR and its partner are known to bind. Disobind cannot reliably distinguish binders from non-binders (see Supplementary text). Third, the accuracy of the predictions depends on the ability of the pLM to provide accurate representations of the IDR and its partner. Finally, our method cannot be used to assess the effects of post-translational modifications as the pLM used does not distinguish post-translationally modified amino acids.


### Challenges

  Predicting contact maps and interface residues for IDRs in a complex with a partner is challenging. First, a limited number of experimental structures are available for IDRs in complexes (Jahn et al., 2024; J. Zhang et al., 2024). Second, IDRs adopt an ensemble of conformations, and the available structures may only partially capture this conformational diversity. Third, inter-protein contact maps are typically sparse, with only a few residues forming contacts. Although not sufficient, we gather all available structures for IDRs in complexes. With these, we create a dataset of merged binary complexes for training Disobind. Using merged binary complexes helps overcome the issue of sparsity associated with training on an ensemble of contact maps. Predicting interface residues instead of contacts and using coarse-graining helps overcome the challenges associated with the multivalency of IDRs and the sparsity of inter-protein contact maps (Fig 3a).


### Performance of AF2 and AF3

  Several studies indicate that the low-confidence regions in AF2 predictions overlap with the presence of disordered regions, although these low-confidence regions cannot be considered as a representative conformation for IDRs (Alderson et al., 2022; Escobedo et al., 2023; Ruff & Pappu, 2021). Despite the success of AF2 and AF3 for ordered proteins, their performance on IDR complexes remains limited, as shown in benchmarks (Alderson et al., 2022; Bret et al., 2024; C. Y. Lee et al., 2024; Omidi et al., 2024). This could potentially be due to several reasons. First, IDRs show low sequence conservation, resulting in poor quality MSAs (Holehouse & Kragelund, 2024). Second, IDRs are best described by an ensemble of conformations which capture their dynamic behaviour. AF2 and AF3, however, predict a single best structure, which may not be representative of IDRs (Ruff & Pappu, 2021). Further, the diversity of binding modes of IDRs in complexes, e.g., DOR and DDR, poses a significant challenge. In particular, AlphaFold predictions are less accurate for fuzzy complexes or DDR complexes, compared to DOR complexes (Alderson et al., 2022; Omidi et al., 2024). Third, the confidence metrics of AF2 and AF3 may not reliably reflect structural accuracy (Elofsson, 2025; Guan & Keating, 2025). Particularly, the presence of IDRs is known to negatively impact the confidence metrics (Dunbrack, 2025; Magana & Kovalevskiy, 2024; Omidi et al., 2024; Varga et al., 2024).

  On our OOD test set, most of the predictions from AF2 and AF3 had low confidence, with an ipTM score lower than 0.75. Notably, compared to AF2, AF3 predictions were less accurate and had lower confidence. It is possible that the confidence cutoffs used for ordered proteins do not apply to IDRs and predictions involving the latter might require different cutoffs. A recent benchmark suggested a lower ipTM score of 0.4 for assessing interfaces involving IDRs (Omidi et al., 2024). Further, the per-residue metrics such as pLDDT and PAE may be more relevant for assessing AF2 and AF3 predictions than the global metrics such as the ipTM and pTM, as is also suggested by a recent benchmark and the AF3 guide (Abramson et al., 2024; Omidi et al., 2024). Thus, AlphaFold confidence metrics might need to be interpreted differently for IDRs, and this requires rigorous benchmarking (Chakravarty et al., 2025; Kim et al., 2024; Varga et al., 2024).

  AF2 has been shown to successfully predict the structures of protein-peptide complexes and domain-motif complexes, although the results could be sensitive to the sequence fragment size, fragment delimitation, and the MSA alignment mode (Bret et al., 2024; B.-G. Lee et al., 2020). Some studies show that MSA subsampling (del Alamo et al., 2022), clustering the sequences in the MSA (Wayment-Steele et al., 2024), and using sliding fragments as input to AF2 (Bret et al., 2024) result in better predictions and can be used to predict multiple conformations. In our comparison, we did not explore these strategies, though they may improve the model predictions.


### Future Directions

  One of the major roadblocks in training methods such as Disobind is the lack of data. More experimental structures of IDRs in complexes would be valuable (Jahn et al., 2024; J. Zhang et al., 2024). Whereas IDR ensembles derived from MD can be used, generating MD ensembles for IDR complexes is computationally expensive and challenging, and the existing databases like PED contain very few such ensembles for IDR complexes (Ghafouri et al., 2024; Majila et al., 2024). Alternatively, deep generative models can be used to generate ensembles for IDRs in complexes (Janson & Feig, 2024; Majila et al., 2024; Mansoor et al., 2024). However, the current methods are limited to generating ensembles for monomers.

  Disobind and similar methods can be further enhanced by improving the existing pLMs to provide better representations for IDRs plausibly by incorporating physical priors and/or structural information (Majila et al., 2024; Rogers et al., 2023; Wang et al., 2024). Additionally, these methods could be extended to predict protein-protein interactions (PPIs) involving IDRs, and aid in the design of IDR binders. The behaviour of IDRs within cells remains largely unexplored. Methods like Disobind are expected to facilitate an improved understanding of the interactions, function, and modulation of IDRs.


<!-- PAPER_BREAK -->

---

<a id="pmid-40286477-amino-acid-sequence-based-idr-classification-using-ensemble-machine-learning-and-quantum-neural-networks"></a>

# PMID 40286477 - Amino acid sequence-based IDR classification using ensemble machine learning and quantum neural networks.

## Title

Amino acid sequence-based IDR classification using ensemble machine learning and quantum neural networks.

## Pmid

40286477

## Keywords

Deep neural network, Intrinsically disordered region, Machine learning, Quantum neural network

## Abstract

Biologically traditional methods, such as the Uversky plot, which rely on hydrophobicity and net charge, have inherent limitations in accurately distinguishing intrinsically disordered regions (IDRs) from ordered protein regions. To overcome these constraints, we propose a novel ensemble framework integrating Machine Learning (ML), Deep Neural Networks (DNN), and Quantum Neural Networks (QNN) to enhance IDR classification accuracy. Notably, this study is the first to employ QNNs for IDR classification, leveraging quantum entanglement to model intricate feature interactions. Amino acid sequences were analyzed to extract biophysical features, including charge distribution, hydrophobicity, and structural properties, which served as inputs for the predictive models. ML was utilized for independent feature learning, DNN for hierarchical interaction modeling, and QNN for capturing high-order dependencies. Our meta-model demonstrated an accuracy of 0.85, surpassing individual classifiers and highlighting the importance of buried amino acids and feature interactions between scaled hydrophobicity and large, buried, and charged residues. This study advances computational protein science by demonstrating the applicability of QNNs in bioinformatics and establishing a robust framework for IDR classification.


<!-- PAPER_BREAK -->

---

<a id="pmid-41378882-impact-of-intrinsically-disordered-regions-and-functional-disorder-hotspots-in-the-human-kinome"></a>

# PMID 41378882 - Impact of intrinsically disordered regions and functional disorder hotspots in the human kinome.

## Title

Impact of intrinsically disordered regions and functional disorder hotspots in the human kinome.

## Pmid

41378882

## Keywords

IDR conformation map, functional disorder hotspots, human kinome, long short-term memory

## Abstract

A ubiquitous and reversible phosphorylation is important for molecular signaling cascades, regulated by the transient interaction of protein kinases. The coupled folding and phosphorylation determining substrate specificity re-calibrates the interactive environment of intrinsically disordered regions (IDRs). There are over 50 computational methods for predicting IDRs in the proteome, yet achieving an accurate depiction remains an ongoing challenge. In this study, we present a standardized and kinase-centric approach for IDR prediction within the human kinome, employing a long short-term memory deep learning framework that achieves a high predictive performance (AUC = 0.97). The web server is now publicly accessible at: https://ciods.in/kindisorder. Our workflow begins with proteome-wide IDR prediction and proceeds with the categorization of short and long IDR segments, followed by an in-depth analysis of their distribution relative to the kinase domain regulatory core. We evaluated the conservation of these IDRs across all 137 human kinase families, computing a trend-setting conservation index to identify both conserved and variable disorder patterns. Through this framework, we uncovered 1039 functional disorder region hotspots that correlate with dynamic conformational shifts, phosphorylation sites, functional motif enrichment, and mutation impact embedded within IDRs. To further validate their regulatory significance, we conducted biophysical profiling of conserved and variable IDRs. Finally, we developed a structural integrity framework to link these IDRs to their influence on intrinsic signaling cascades and substrate specificity. This study offers a comprehensive functional characterization of IDRs in the human kinome, providing a valuable resource for exploring kinase regulation and opportunities in drug repurposing.

## Conclusion

This study presents a high-confidence prediction of IDR in the human kinome. Current prediction tools heavily compress the proteome-level insights. They also struggle with accurately identifying short disorder regions, which play crucial roles in molecular recognition [60]. Our model outperforms flDPnn, flDPnn2, ADOPT, DisoFLAG, AIUPred in AUC, and ROC metrics, highlighting its predictive strength. While flDPnn2 surpasses many existing tools, most models struggle with long sequences, introducing noise at the proteome scale. Many cannot handle proteins over 5000 amino acids, limiting their use in large kinases like OBSCN and TTK. In contrast, our approach provides a scalable, accurate solution for disorder prediction across diverse kinases.

To further probe the functional significance of kinase IDRs, we categorized the IDRs in domain-specific categorization involving a logical analysis bound to short and long stretch IDRs. It is imperative that molecular recognition property identification and interpretation of flexibility segments in understudied kinases. Using a combination of functional motifs, distance metrics, mutation probability, MultiPTM context, and reported phosphosites evaluated how functional DR hotspots alter kinase modularity. This allowed us to explore how sequence variations within IDRs contribute to disease-associated dysfunctions or regulatory modifications.

Statistical comparisons between short and long disorder scores revealed distinct distribution patterns, suggesting that long disordered stretches are more conserved and potentially more functionally relevant. Detailed mapping across kinase domains further highlighted that long disorder regions often lie outside the canonical kinase domains, frequently overlapping with hotspots of post-translational modifications and functional motifs, thereby implying regulatory roles. Conversely, short disorder regions, although present both within and outside domains, showed less pronounced conservation. Finally, mutation analysis of a representative AKT kinase demonstrated the structural and energetic impact of a point mutation within a disordered region, with a moderate shift in ΔG and significant transition probability, underscoring the potential of disorder-associated residues to influence protein dynamics. Together, these findings underscore the critical regulatory and functional roles of IDRs in kinase biology and provide a valuable framework for further exploration of disorder-mediated signaling and disease-associated mutations. While our free energy-based model effectively explains the impact of many site-specific mutations by integrating sequence, secondary structure, and side chain information, it is not without limitations. In particular, the model may underperform for mutations involving long-range allosteric effects, regions with sparse structural or evolutionary data, or those influenced by post-translational modifications. These cases highlight areas for future refinement and the potential integration of additional contextual or structural dynamics data.


<!-- PAPER_BREAK -->

---

<a id="pmid-41534519-disobind-a-sequence-based-partner-dependent-contact-map-and-interface-residue-predictor-for-intrinsically-disordered-regions"></a>

# PMID 41534519 - Disobind: A sequence-based, partner-dependent contact map and interface residue predictor for intrinsically disordered regions.

## Title

Disobind: A sequence-based, partner-dependent contact map and interface residue predictor for intrinsically disordered regions.

## Pmid

41534519

## Keywords

DL, IDP, IDR, deep learning, intrinsically disordered proteins, intrinsically disordered regions, pLMs, protein language model, protein structure

## Abstract

Intrinsically disordered proteins or regions (IDPs or IDRs) adopt diverse binding modes with different partners, ranging from coupled folding and binding to fuzzy binding and fully disordered binding. Characterizing IDR interfaces is challenging both experimentally and computationally. State-of-the-art tools such as AlphaFold multimer and AlphaFold3 can be used to predict IDR binding sites, although they are less accurate at their benchmarked confidence cutoffs. Here, we developed Disobind, a deep-learning method that predicts inter-protein contact maps and interface residues for an IDR and its partner, given their sequences. It uses sequence embeddings from the ProtT5 protein language model. Disobind outperforms state-of-the-art interface predictors for IDRs. It also outperforms AlphaFold multimer and AlphaFold3 at multiple confidence cutoffs. Combining Disobind and AlphaFold-multimer predictions further improves performance. In contrast to current methods, Disobind considers the context of the binding partner and does not depend on structures and multiple sequence alignments. Its predictions can be used to localize IDRs in large assemblies and characterize IDR-mediated interactions.

