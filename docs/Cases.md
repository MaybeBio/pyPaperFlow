# Test Cases

[Back to README](../README.md)

---

## Index

[Case1: Get PMIDs from Query](#-case-1-get-pmids-from-query)

[Case2: Fetch Metadata for pubmed papers from query or PMIDs list](#-case-2-fetch-metadata-for-pubmed-papers-from-query-or-pmids-list)

[Case3: Fetch full text data for pubmed papers from PMIDs list](#-case-3-fetch-full-text-data-for-pubmed-papers-from-pmids-list)

[Case4: Fetch full paper data (including metadata and full text data) for pubmed papers from PMIDs list](#-case-4-fetch-full-paper-data-including-metadata-and-full-text-data-for-pubmed-papers-from-pmids-list)

[Case5: Prepare batch Markdown-formatted paper data for downstream LLMs.](#-case-5-prepare-batch-markdown-formatted-paper-data-for-downstream-llms)



## 🧬 Case 1: Get PMIDs from Query

run the command:
```bash
paperflow pubmed-search "alphafold3 AND conformation AND ensemble" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test
```
the log shows:
```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
Now searching PubMed with query [alphafold3 AND conformation AND ensemble] at [2026-05-02 17:03:20] ...
found 19 related articles about [alphafold3 AND conformation AND ensemble] at [2026-05-02 17:03:22] ...
Retrieving 19 PMIDs from history server at [2026-05-02 17:03:22] ...
Fetching PMIDs 1 to 19 at [2026-05-02 17:03:22] ...
  -> Retrieved 19 PMIDs in this batch.
Total PMIDs retrieved: 19 out of 19 at [2026-05-02 17:03:23] ...
Found 19 PMIDs.
['41914502', '41779774', '41639320', '41502950', '41478913', '41432299', '41249430', '41147497', '41047853', '41014267', '40950168', '40938899', '40714407', '40549150', '40490178', '39574676', '39186607', '38996889', '38995731']
PMIDs saved to ./test/pubmed_searched_ids_2026-05-02_17-03-23.txt.
```
As you can see, we will print the PMIDs list for you and save it in a text file which can be used further.

⚠️ We also recommend using the `Search & Save plugin` on the PubMed webpage to obtain the PMID list for subsequent use.

![alt text](../figs/pubmed.png)

## 🧬 Case 2: Fetch Metadata for pubmed papers from query or PMIDs list

If you do not have detailed PMID list and want to fetch meta information from query, run the command:
```bash
paperflow pubmed-meta -q "alphafold3 AND conformation AND ensemble" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test/alphafold3_ensemble_meta
```

the log shows:
```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
Fetching papers for query: alphafold3 AND conformation AND ensemble
Now searching PubMed with query [alphafold3 AND conformation AND ensemble] at [2026-05-02 17:07:55] ...
found 19 related articles about [alphafold3 AND conformation AND ensemble] at [2026-05-02 17:07:56] ...
Fetching articles 1 to 19 at [2026-05-02 17:07:56] ...
  -> Retrieved 19 Medline records and 19 Xml articles. Please check whether they equal and the efetch number here with esearch count.
  -> Deep mining 5 types of internal connections for 19 PMIDs at [2026-05-02 17:07:59] ...
     Fetching pubmed_pubmed_refs from pubmed for 19 PMIDs at [2026-05-02 17:07:59] ...
     Fetching pubmed_pubmed from pubmed for 19 PMIDs at [2026-05-02 17:08:02] ...
     Fetching pubmed_pubmed_reviews from pubmed for 19 PMIDs at [2026-05-02 17:08:05] ...
     Fetching pubmed_pmc from pmc for 19 PMIDs at [2026-05-02 17:08:08] ...
     Fetching pubmed_pubmed_citedin from pubmed for 19 PMIDs at [2026-05-02 17:08:12] ...
  -> Fetching external LinkOuts (Datasets, Full Text, etc.) for 19 PMIDs at [2026-05-02 17:08:15] ...
  -> Saved 41914502 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41914502/41914502_meta.json
  -> Saved 41779774 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41779774/41779774_meta.json
  -> Saved 41639320 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41639320/41639320_meta.json
  -> Saved 41502950 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41502950/41502950_meta.json
  -> Saved 41478913 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41478913/41478913_meta.json
  -> Saved 41432299 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41432299/41432299_meta.json
  -> Saved 41249430 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/41249430/41249430_meta.json
  -> Saved 41147497 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41147497/41147497_meta.json
  -> Saved 41047853 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41047853/41047853_meta.json
  -> Saved 41014267 metadata to ./test/alphafold3_ensemble_meta/pubmed/2026/41014267/41014267_meta.json
  -> Saved 40950168 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/40950168/40950168_meta.json
  -> Saved 40938899 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/40938899/40938899_meta.json
  -> Saved 40714407 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/40714407/40714407_meta.json
  -> Saved 40549150 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/40549150/40549150_meta.json
  -> Saved 40490178 metadata to ./test/alphafold3_ensemble_meta/pubmed/2025/40490178/40490178_meta.json
  -> Saved 39574676 metadata to ./test/alphafold3_ensemble_meta/pubmed/2024/39574676/39574676_meta.json
  -> Saved 39186607 metadata to ./test/alphafold3_ensemble_meta/pubmed/2024/39186607/39186607_meta.json
  -> Saved 38996889 metadata to ./test/alphafold3_ensemble_meta/pubmed/2024/38996889/38996889_meta.json
  -> Saved 38995731 metadata to ./test/alphafold3_ensemble_meta/pubmed/2024/38995731/38995731_meta.json
```

you can check the result here: [alphafold3_ensemble_meta](../test/alphafold3_ensemble_meta/)

As shown above, a `/pubmed` subfolder will be automatically created under your output directory, with all metadata JSON files saved inside this folder.

Otherwise, if you have detailed PMID list,
run the command below:
```bash
# here we use search list in case1 as an example
paperflow pubmed-meta -f ./test/pubmed_searched_ids_2026-05-02_17-03-23.txt  --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test/alphafold3_ensemble_meta_try2
```

the log shows the same way:

```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
Fetching 19 papers from file /data2/pyPaperFlow/test/pubmed_searched_ids_2026-05-02_17-03-23.txt.
Total PMIDs to fetch: 19 at [2026-05-02 17:14:40] ...
Fetching articles 1 to 19 (PMID: ['41914502', '41779774', '41639320', '41502950', '41478913', '41432299', '41249430', '41147497', '41047853', '41014267', '40950168', '40938899', '40714407', '40549150', '40490178', '39574676', '39186607', '38996889', '38995731']) at [2026-05-02 17:14:40] ...
  -> Retrieved 19 Medline records and 19 Xml articles. Please check whether they equal and whether they match the number of this batch.
  -> Deep mining 5 types of internal connections for 19 PMIDs at [2026-05-02 17:14:43] ...
     Fetching pubmed_pmc from pmc for 19 PMIDs at [2026-05-02 17:14:43] ...
     Fetching pubmed_pubmed_citedin from pubmed for 19 PMIDs at [2026-05-02 17:14:46] ...
     Fetching pubmed_pubmed_refs from pubmed for 19 PMIDs at [2026-05-02 17:14:49] ...
     Fetching pubmed_pubmed_reviews from pubmed for 19 PMIDs at [2026-05-02 17:14:53] ...
     Fetching pubmed_pubmed from pubmed for 19 PMIDs at [2026-05-02 17:14:57] ...
  -> Fetching external LinkOuts (Datasets, Full Text, etc.) for 19 PMIDs at [2026-05-02 17:15:01] ...
  -> Saved 41914502 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41914502/41914502_meta.json
  -> Saved 41779774 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41779774/41779774_meta.json
  -> Saved 41639320 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41639320/41639320_meta.json
  -> Saved 41502950 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41502950/41502950_meta.json
  -> Saved 41478913 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41478913/41478913_meta.json
  -> Saved 41432299 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41432299/41432299_meta.json
  -> Saved 41249430 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/41249430/41249430_meta.json
  -> Saved 41147497 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41147497/41147497_meta.json
  -> Saved 41047853 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41047853/41047853_meta.json
  -> Saved 41014267 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2026/41014267/41014267_meta.json
  -> Saved 40950168 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/40950168/40950168_meta.json
  -> Saved 40938899 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/40938899/40938899_meta.json
  -> Saved 40714407 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/40714407/40714407_meta.json
  -> Saved 40549150 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/40549150/40549150_meta.json
  -> Saved 40490178 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2025/40490178/40490178_meta.json
  -> Saved 39574676 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2024/39574676/39574676_meta.json
  -> Saved 39186607 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2024/39186607/39186607_meta.json
  -> Saved 38996889 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2024/38996889/38996889_meta.json
  -> Saved 38995731 metadata to ./test/alphafold3_ensemble_meta_try2/pubmed/2024/38995731/38995731_meta.json
```

We store the meta data of the paper in a json file.
 
One example [PMID 41249430](./test/alphafold3_ensemble_meta/pubmed/2025/41249430/41249430_meta.json) listed as below:

```python
{
    "content": {
        "abstract": "AlphaFold2 and AlphaFold3 have revolutionized protein structure prediction by enabling high-accuracy structure predictions for most single-chain proteins. However, obtaining high-quality predictions for difficult targets with shallow or noisy multiple sequence alignments and complicated multi-domain architectures remains challenging. We present MULTICOM4, an integrative structure prediction system that uses diverse MSA generation, large-scale model sampling, and an ensemble model quality assessment strategy to improve model generation and ranking of AlphaFold2 and AlphaFold3. In the 16th Critical Assessment of Techniques for Protein Structure Prediction, our predictors built on MULTICOM4 ranked among the top out of 120 predictors in tertiary structure prediction and outperformed a standard AlphaFold3 predictor. Our best predictor achieved an average TM-score of 0.902 for 84 CASP16 domains, with top-1 predictions reaching high accuracy (TM-score>0.9) for 73.8% and correct folds (TM-score>0.5) for 97.6% of domains. For best-of-top-5 predictions, all domains were correctly folded. The results show that MSA engineering using different sequence databases, alignment tools, and domain segmentation along with extensive model sampling is critical to generate accurate structural models. Combining complementary QA methods with model clustering further improves ranking reliability. These advances provide practical strategies for modeling difficult single-chain proteins in structural biology and drug discovery.",
        "keywords": [],
        "mesh_terms": [
            "*Computational Biology/methods",
            "Protein Folding",
            "Models, Molecular",
            "*Protein Structure, Tertiary",
            "*Proteins/chemistry",
            "Sequence Alignment/methods",
            "*Software",
            "Sequence Analysis, Protein/methods",
            "Algorithms"
        ],
        "pub_types": [
            "Journal Article"
        ]
    },
    "contributors": {
        "medline": {
            "affiliations": [
                "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA.",
                "NextGen Precision Health, University of Missouri, Columbia, MO, USA.",
                "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA.",
                "NextGen Precision Health, University of Missouri, Columbia, MO, USA.",
                "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA. chengji@missouri.edu.",
                "NextGen Precision Health, University of Missouri, Columbia, MO, USA. chengji@missouri.edu."
            ],
            "auids": [
                "ORCID: 0000-0003-0305-2853"
            ],
            "full_names": [
                "Liu, Jian",
                "Neupane, Pawan",
                "Cheng, Jianlin"
            ],
            "short_names": [
                "Liu J",
                "Neupane P",
                "Cheng J"
            ]
        },
        "xml": [
            {
                "affiliations": [
                    "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA.",
                    "NextGen Precision Health, University of Missouri, Columbia, MO, USA."
                ],
                "full_name": "Liu, Jian",
                "identifiers": [],
                "short_name": "Liu J"
            },
            {
                "affiliations": [
                    "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA.",
                    "NextGen Precision Health, University of Missouri, Columbia, MO, USA."
                ],
                "full_name": "Neupane, Pawan",
                "identifiers": [],
                "short_name": "Neupane P"
            },
            {
                "affiliations": [
                    "Department of Electrical Engineering & Computer Science, University of Missouri, Columbia, MO, USA. chengji@missouri.edu.",
                    "NextGen Precision Health, University of Missouri, Columbia, MO, USA. chengji@missouri.edu."
                ],
                "full_name": "Cheng, Jianlin",
                "identifiers": [
                    "0000-0003-0305-2853"
                ],
                "short_name": "Cheng J"
            }
        ]
    },
    "identity": {
        "doi": "10.1038/s42003-025-08960-6",
        "pmid": "41249430",
        "title": "Boosting AlphaFold protein tertiary structure prediction through MSA engineering and extensive model sampling and ranking in CASP16."
    },
    "links": {
        "cites": [
            "41178755",
            "40672254"
        ],
        "entrez": {},
        "external": [
            {
                "attribute": "free resource",
                "category": "Full Text Sources",
                "linkname": "",
                "provider": "Nature Publishing Group",
                "url": "https://doi.org/10.1038/s42003-025-08960-6"
            },
            {
                "attribute": "free resource",
                "category": "Full Text Sources",
                "linkname": "",
                "provider": "PubMed Central",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/pmid/41249430/"
            },
            {
                "attribute": "free resource",
                "category": "Research Materials",
                "linkname": "",
                "provider": "NCI CPTC Antibody Characterization Program",
                "url": "https://antibodies.cancer.gov/detail/CPTC-TOP1-1"
            }
        ],
        "pmc": [
            "12623963"
        ],
        "refs": [
            "40799498",
            "40452318",
            "39123049",
            "38718835",
            "38167654",
            "37949999",
            "37679431",
            "36927031",
            "36734597",
            "34873061",
            "34453465",
            "34291486",
            "34282049",
            "34265844",
            "31942072",
            "31696235",
            "31676016",
            "31399549",
            "31235882",
            "30395287",
            "29959318",
            "29228193",
            "29228185",
            "27899574",
            "25391399",
            "24225321",
            "23047561",
            "22198341",
            "20718988",
            "18542861",
            "11159328"
        ],
        "review": [
            "41249430",
            "38316555",
            "38986287",
            "40973394",
            "39970826",
            "40332289"
        ],
        "similar": [
            "41249430",
            "40661500",
            "40585263",
            "40452318",
            "40161604",
            "41170922",
            "41014267",
            "40820259",
            "40851426",
            "40501681",
            "40762404",
            "41147497",
            "40751131",
            "37650367",
            "19077267",
            "40847537",
            "17553833",
            "40874652",
            "40799498",
            "41104652",
            "34599769",
            "37949999",
            "34382712",
            "26369671",
            "40502139",
            "38316555",
            "34331351",
            "31344267",
            "19701941",
            "41313605",
            "20066664",
            "19777061",
            "34240477",
            "34162922",
            "30985027",
            "28093407",
            "38986287",
            "24637808",
            "41165252",
            "40950168",
            "14579329",
            "34455641",
            "37293073",
            "37679431",
            "40195868",
            "19722267",
            "40810260",
            "40488225",
            "25431331",
            "28748648",
            "41047853",
            "37565699",
            "18452616",
            "34291486",
            "18487301",
            "16187361",
            "26445311",
            "41201924",
            "16187348",
            "18215316",
            "37321965",
            "41257887",
            "24018415",
            "34884640",
            "41081541",
            "35034173",
            "39052676",
            "29082551",
            "17452345",
            "22069035",
            "41454828",
            "41325379",
            "40778521",
            "31365149",
            "31471916",
            "40973394",
            "39970826",
            "15359422",
            "27028541",
            "17570145",
            "14579328",
            "29139163",
            "22168237",
            "40332289",
            "21301031",
            "23812990",
            "40696837",
            "33850214",
            "26713437",
            "41045049",
            "26343917",
            "38913900",
            "31918654",
            "40067116",
            "20470364",
            "15939584",
            "22545707",
            "17680686",
            "41261173",
            "31634369"
        ],
        "text_mined": []
    },
    "metadata": {
        "entrez_date": "2025/11/18 00:28",
        "fetched_at": "2026-05-02 15:17:48"
    },
    "source": {
        "journal_abbrev": [
            "Commun Biol"
        ],
        "journal_title": [
            "Communications biology"
        ],
        "pub_date": "2025 Nov 17",
        "pub_types": [
            "Journal Article"
        ],
        "pub_year": "2025"
    }
}

```

![alt text](./figs/41249430.png)


## 🧬 Case 3: Fetch full text data for pubmed papers from PMIDs list

If you only need to fetch the full text from PMIDs — where the `full text refers to the main body of a paper` (the complete textual content equivalent to that parsed from PDF files) — you can simply run
```bash
# we choose pmid 39570595 here as an example
paperflow pubmed-content -p 39570595  --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test/full_text
```

the log shows
```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
Downloading full texts for 1 PMIDs from file provided PMIDs.
Fetching full text for 1 Pubmed articles at [2026-05-02 17:25:41] ...
 -> Converting Pubmed articles 1 to 1 (PMID : ['39570595']) to PMC IDs at [2026-05-02 17:25:41] ...
  -> Mapped 1 out of 1 PMIDs to valid PMC IDs. Downloading full text XML for these PMC IDs at [2026-05-02 17:25:42] ...
  -> Saved XML to ./test/full_text/pubmed/2024/39570595/39570595_content.xml
  -> Saved parsed JSON to ./test/full_text/pubmed/2024/39570595/39570595_content.json
  -> Saved parsed text to ./test/full_text/pubmed/2024/39570595/39570595_content.md
```

As you can see, for full-text data, we handle it differently from metadata—while metadata is simply stored in a JSON file named `{PMID}_meta.json`, full-text data is output into three files with distinct formats, each serving a specific purpose:
* **{PMID}_content.xml**: Stores the raw XML content retrieved directly from the response, preserving the original data structure.
* **{PMID}_content.json**: Contains detailed, structured full-text content. This format allows for direct extraction of specific sections (e.g., introduction, results, discussion), making it ideal for quick exploration or targeted analysis of particular parts of the text.
* **{PMID}_content.md**: Saves the full text of the paper in Markdown format. Its clean, human-readable structure makes it well-suited for high-throughput summarization tasks using LLMs/AI tools (such as ChatGPT or other preferred models).

`Core Principle: JSON for coding, Markdown for LLM prompting.`


or you can batch download what you want 
```bash
# we use searched_pmids.txt generated by Case1
paperflow download-fulltext  -f ./test/pubmed_searched_ids_2026-05-02_17-03-23.txt  --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test/alphafold_ensemble_content_try3
```

the log shows 
```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
Downloading full texts for 19 PMIDs from file /data2/pyPaperFlow/test/pubmed_searched_ids_2026-05-02_17-03-23.txt.
Fetching full text for 19 Pubmed articles at [2026-05-02 18:09:34] ...
 -> Converting Pubmed articles 1 to 19 (PMID : ['41914502', '41779774', '41639320', '41502950', '41478913', '41432299', '41249430', '41147497', '41047853', '41014267', '40950168', '40938899', '40714407', '40549150', '40490178', '39574676', '39186607', '38996889', '38995731']) to PMC IDs at [2026-05-02 18:09:34] ...
  -> Mapped 10 out of 19 PMIDs to valid PMC IDs. Downloading full text XML for these PMC IDs at [2026-05-02 18:09:36] ...
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2025/40950168/40950168_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2025/40950168/40950168_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2025/40950168/40950168_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2026/41914502/41914502_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2026/41914502/41914502_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2026/41914502/41914502_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2025/41432299/41432299_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2025/41432299/41432299_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2025/41432299/41432299_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2025/40549150/40549150_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2025/40549150/40549150_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2025/40549150/40549150_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2026/41147497/41147497_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2026/41147497/41147497_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2026/41147497/41147497_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2024/39574676/39574676_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2024/39574676/39574676_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2024/39574676/39574676_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2025/41249430/41249430_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2025/41249430/41249430_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2025/41249430/41249430_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2024/38995731/38995731_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2024/38995731/38995731_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2024/38995731/38995731_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2025/40938899/40938899_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2025/40938899/40938899_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2025/40938899/40938899_content.md
  -> Saved XML to ./test/alphafold_ensemble_content_try3/pubmed/2026/41502950/41502950_content.xml
  -> Saved parsed JSON to ./test/alphafold_ensemble_content_try3/pubmed/2026/41502950/41502950_content.json
  -> Saved parsed text to ./test/alphafold_ensemble_content_try3/pubmed/2026/41502950/41502950_content.md

```
as you can imagine, not all pmids have a validated pmc id, you can try other tools for free full text extraction. 


## 🧬 Case 4: Fetch full paper data (including metadata and full text data) for pubmed papers from PMIDs list

Now if you want to get everything of papers you want, not just metadata or full text but BOTH!

You can simply run 
```bash
# from query 
paperflow pubmed-all --query "IDR AND interaction AND deep learning" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY -o ./test/full_paper_test

# from PMID list, same as above
```

for query subcommand, the log shows 
```bash
✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.
=== Step 1: Fetching Metadata ===
Now searching PubMed with query [IDR AND interaction AND deep learning] at [2026-05-02 18:19:06] ...
found 6 related articles about [IDR AND interaction AND deep learning] at [2026-05-02 18:19:07] ...
Fetching articles 1 to 6 at [2026-05-02 18:19:07] ...
  -> Retrieved 6 Medline records and 6 Xml articles. Please check whether they equal and the efetch number here with esearch count.
  -> Deep mining 5 types of internal connections for 6 PMIDs at [2026-05-02 18:19:10] ...
     Fetching pubmed_pubmed from pubmed for 6 PMIDs at [2026-05-02 18:19:10] ...
     Fetching pubmed_pubmed_reviews from pubmed for 6 PMIDs at [2026-05-02 18:19:12] ...
     Fetching pubmed_pmc from pmc for 6 PMIDs at [2026-05-02 18:19:14] ...
     Fetching pubmed_pubmed_refs from pubmed for 6 PMIDs at [2026-05-02 18:19:15] ...
     Fetching pubmed_pubmed_citedin from pubmed for 6 PMIDs at [2026-05-02 18:19:17] ...
  -> Fetching external LinkOuts (Datasets, Full Text, etc.) for 6 PMIDs at [2026-05-02 18:19:19] ...
  -> Saved 41534519 metadata to ./test/full_paper_test/pubmed/2026/41534519/41534519_meta.json
  -> Saved 41378882 metadata to ./test/full_paper_test/pubmed/2025/41378882/41378882_meta.json
  -> Saved 40286477 metadata to ./test/full_paper_test/pubmed/2025/40286477/40286477_meta.json
  -> Saved 39763873 metadata to ./test/full_paper_test/pubmed/2025/39763873/39763873_meta.json
  -> Saved 38701796 metadata to ./test/full_paper_test/pubmed/2024/38701796/38701796_meta.json
  -> Saved 36851914 metadata to ./test/full_paper_test/pubmed/2023/36851914/36851914_meta.json

=== Step 2: Fetching Full Text ===
Fetching full text for 6 Pubmed articles at [2026-05-02 18:19:22] ...
 -> Converting Pubmed articles 1 to 6 (PMID : ['41534519', '41378882', '40286477', '39763873', '38701796', '36851914']) to PMC IDs at [2026-05-02 18:19:22] ...
  -> Mapped 3 out of 6 PMIDs to valid PMC IDs. Downloading full text XML for these PMC IDs at [2026-05-02 18:19:24] ...
  -> Saved XML to ./test/full_paper_test/pubmed/2023/36851914/36851914_content.xml
  -> Saved parsed JSON to ./test/full_paper_test/pubmed/2023/36851914/36851914_content.json
  -> Saved parsed text to ./test/full_paper_test/pubmed/2023/36851914/36851914_content.md
  -> Saved XML to ./test/full_paper_test/pubmed/2025/41378882/41378882_content.xml
  -> Saved parsed JSON to ./test/full_paper_test/pubmed/2025/41378882/41378882_content.json
  -> Saved parsed text to ./test/full_paper_test/pubmed/2025/41378882/41378882_content.md
  -> Saved XML to ./test/full_paper_test/pubmed/2025/39763873/39763873_content.xml
  -> Saved parsed JSON to ./test/full_paper_test/pubmed/2025/39763873/39763873_content.json
  -> Saved parsed text to ./test/full_paper_test/pubmed/2025/39763873/39763873_content.md

=== Step 3: Processing and Saving Metadata ===
  -> Saved 41534519 metadata to ./test/full_paper_test/pubmed/2026/41534519/41534519_meta.json
  -> Extracted 2 URLs from full text for PMID 41378882
  -> Saved 41378882 metadata to ./test/full_paper_test/pubmed/2025/41378882/41378882_meta.json
  -> Saved 40286477 metadata to ./test/full_paper_test/pubmed/2025/40286477/40286477_meta.json
  -> Extracted 2 URLs from full text for PMID 39763873
  -> Saved 39763873 metadata to ./test/full_paper_test/pubmed/2025/39763873/39763873_meta.json
  -> Saved 38701796 metadata to ./test/full_paper_test/pubmed/2024/38701796/38701796_meta.json
  -> Extracted 29 URLs from full text for PMID 36851914
  -> Saved 36851914 metadata to ./test/full_paper_test/pubmed/2023/36851914/36851914_meta.json

```

As shown above, two types of files will be generated: `{PMID}_meta.*` and `{PMID}_content.*`.


## 🧬 Case 5: Prepare batch Markdown-formatted Pubmed paper data for downstream LLMs. 

Once you have retrieved all relevant papers(meta+content) on a specific topic or theme, the next step is to aggregate them into a unified collection.
In this step, we merge all papers with complete metadata and full content into a `paper-level JSON file` for consolidated summarization.
You may also extract designated sections of these papers—such as the abstract, discussion, and conclusion—and compile them into a well-structured Markdown file, which is fully compatible with `downstream LLM text-based parsing tasks`.

```bash
# paper directory from Case4 result
paperflow pubmed-merge-json -i /data2/pyPaperFlow/test/full_paper_test -o /data2/pyPaperFlow/test/full_paper_test --jsonl -s /data2/pyPaperFlow/test/full_paper_test
```
the log shows

```
✅ Please check the merged pubmed JSON/JSONL file at /data2/pyPaperFlow/test/full_paper_test and the merge statistics file at /data2/pyPaperFlow/test/full_paper_test. Also, a JSON file per paper is created within the PMID subfolders.
```

You can access the merged JSONL file [here](./test/full_paper_test/full_paper_test_2026-05-05_22-02-51.jsonl), where each line corresponds to one paper in JSON format. The statistical results are also available [here](./test/full_paper_test/full_paper_test_stats_2026-05-05_22-02-51.json).

In statistical JSON file, you can see PMID `"38701796","40286477","41534519"` paper is content-missing. 

For these papers, we provide a DOI-based PDF retrieval module, along with another module that parses PDF files into Markdown format, which is fully compatible with the aforementioned `{PMID}.json` files.

Additionally, each paper has a corresponding `{PMID}.json` file containing both metadata and full content information. 

Next, you can use the merged JSONL/JSON file to extract specific sections of interest (e.g., abstract, discussion, conclusion) and compile them into a well-structured Markdown file for downstream LLM text-based parsing tasks.

But you need to provide a configuration file to specify the sections to extract.

Here is an example configuration file:
```bash
metadata_fields:
  - identity.title
  - identity.pmid
  - identity.doi
  - content.keywords
  - content.mesh_terms
  - content.pub_types
  - content.abstract # abstract in metadata first, fall back in content sections(deprecated)
  - contributors.medline
  - contributors.xml
  - links.cites
  - links.entrez
  - links.external
  - links.pmc
  - links.refs
  - links.review
  - links.similar
  - links.text_mined
  - metadata.entrez_date
  - metadata.fetched_at
  - source.journal_abbrev
  - source.journal_title
  - source.pub_date
  - source.pub_types
  - source.pub_year

content_sections:
  - abstract
  - introduction
  - methods
  - results
  - discussion
  - conclusion
  - supplementary
  - availability
  - funding
  - acknowledgements
  - author_contributions

```

As you can see, the `metadata_fields` are actually the keys in the metadata JSON file, and you can specify which fields to extract based on your needs. 
However, the `content_sections` differ from one another, as we extract hierarchical nodes from XML files.
We initially parse distinct nodes in XML format. Fortunately, the `first-level nodes` under the body hierarchy closely correspond to `standard academic sections including abstract, introduction, methods, results, discussion, and conclusion`. For this reason, we directly designate these nodes as content_sections for selective content extraction.
With `simple regular expression matching and string manipulation`, we can realize the mapping and extraction of unified content_sections across different texts.
You may customize which parts to extract by configuring the content_sections field in the `config.yaml` file.
If no custom YAML configuration file is provided, all the above content_sections will be extracted by default, i.e., the full text of the main content.


You can try the following command:
```bash
paperflow pubmed-extract-md -i ./test/full_paper_test/full_paper_test_2026-05-05_22-02-51.jsonl -o ./test/full_paper_test/test_no_yaml.md 

```

the log shows
```bash
Successfully exported 6 papers to /data2/pyPaperFlow/test/full_paper_test/test_no_yaml.md
```

the output markdown file is [here](./test/full_paper_test/test_no_yaml.md)

Each paper is separated by `---` in the markdown file, which is suitable for downstream LLM text-based parsing tasks.

If you want to use a custom YAML configuration file(an example is in [test.yaml](./test/full_paper_test/test.yaml)), you can run the following command:
```bash
paperflow pubmed-extract-md -i ./test/full_test_20_test_2026-05-05_22-02-51.jsonl -o ./test/full_paper_test/test_with_yaml.md -c ./test.yaml
```

the output markdown file is [here](./test/full_paper_test/test_with_yaml.md), where only the sections specified in the YAML file are extracted and compiled into the markdown file.

And next, you can use the generated markdown file for downstream LLM-based summarization or other text-based parsing tasks (e.g. Summarize over the conclusion and discussion sections of all these papers to put forward a research question).


## 🧬 Case 6: Fetch full text data for unaccessible papers based on DOI - fetch PDF then parsing it into markdown 

For those papers that are content-missing (Pubmed papers missing PMC links, or papers with DOI in other databases but no PDF available), we provide a DOI-based PDF retrieval module, along with another module that parses PDF files into Markdown format.

```python 
❯ paperflow paper-fetch --help
usage: paper-fetch [-h] [--title TITLE] [--batch FILE] [--out DIR] [--dry-run]
                   [--format {json,text}] [--pretty] [--stream] [--overwrite]
                   [--idempotency-key KEY] [--timeout SECONDS] [--version]
                   [doi]

Fetch legal open-access PDFs by DOI via Unpaywall, Semantic Scholar, arXiv, PMC, and bioRxiv/medRxiv.

positional arguments:
  doi                   DOI to fetch (e.g. 10.1038/s41586-020-2649-2). Use '-' to
                        read from stdin.

options:
  -h, --help            show this help message and exit
  --title TITLE         paper title; resolved to a DOI via Crossref before download.
                        Mutually exclusive with positional DOI / --batch.
  --batch FILE          file with one DOI per line for bulk download. Use '-' to read
                        from stdin.
  --out DIR             output directory (default: pdfs)
  --dry-run             resolve sources without downloading; preview the PDF URL and
                        filename
  --format {json,text}  output format. json for agents, text for humans. Default:
                        json when stdout is not a TTY, text otherwise.
  --pretty              pretty-print JSON output (2-space indent)
  --stream              emit one NDJSON result per line on stdout as each DOI
                        resolves (batch mode)
  --overwrite           re-download even if the destination file already exists
  --idempotency-key KEY
                        safe-retry key; re-running with the same key replays the
                        original envelope from <out>/.paper-fetch-idem/
  --timeout SECONDS     HTTP timeout in seconds per request (default: 30)
  --version             show program's version number and exit

exit codes:
  0  all DOIs resolved successfully
  1  unresolved (some DOIs had no OA copy; no transport failure)
  3  validation error (bad arguments)
  4  transport error (network / download / IO failure; retryable class)

subcommands:
  schema                 print the machine-readable CLI schema and exit (no network)

stdin:
  paper-fetch -          read a single DOI from stdin
  paper-fetch --batch -  read DOIs line-by-line from stdin

output:
  stdout emits one JSON object per invocation (NDJSON with --stream).
  stderr emits NDJSON progress events when --format json, prose when --format text.
  stdout format auto-detects TTY: json when piped/captured, text in a terminal.

examples:
  paper-fetch 10.1038/s41586-020-2649-2
  paper-fetch 10.1038/s41586-020-2649-2 --dry-run
  paper-fetch --batch dois.txt --out ./papers --format text
  echo 10.1038/s41586-020-2649-2 | paper-fetch --batch -
  paper-fetch schema
```

We recommend you using only arguments below for pdf fetching and leave the rest arguments as default values.
```
--title
--batch
--out
--dry-run
--timeout
```

Here we use paper [IDPFold](https://pubmed.ncbi.nlm.nih.gov/41082321/) as an example, its DOI is `10.1002/advs.202511636`

so you can run the following command to fetch its PDF file:
```bash
paperflow paper-fetch  --out ./test/Other_database --timeout 30  10.1002/advs.202511636
```

the log shows
```bash
==> 10.1002/advs.202511636
  [unpaywall] trying…
  [unpaywall] no PDF
  [semantic_scholar] trying…
  [semantic_scholar] no PDF
  [europe_pmc] https://europepmc.org/articles/PMC12752595?pdf=render
  [pmc] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12752595/pdf/
  saved → /data2/pyPaperFlow/test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.pdf
[europe_pmc] 10.1002/advs.202511636 → /data2/pyPaperFlow/test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.pdf  (saved)

1/1 succeeded  (0 failed)
```

the output pdf file is [here](./test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.pdf)

After you get the PDF file, you can use the following command to parse it into Markdown format:
```bash
paperflow pdf2md 

```

a typical log shows the following:
```
2026-05-09 09:12:55.360 | INFO     | mineru.cli.client:run_orchestrated_cli:874 - Started local mineru-api at http://127.0.0.1:48227
2026-05-09 09:12:56.342 | INFO     | __main__:create_app:260 - Request concurrency limited to 3
Start MinerU FastAPI Service: http://127.0.0.1:48227
API documentation: http://127.0.0.1:48227/docs
INFO:     Started server process [2957300]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://127.0.0.1:48227 (Press CTRL+C to quit)
2026-05-09 09:12:56.364 | INFO     | mineru.cli.client:run_planned_task:771 - Submitting batch 1/1 | 1 document, 12 pages in this batch | 12 pages total | task#1 [Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En]
2026-05-09 09:12:57.588 | INFO     | mineru.backend.pipeline.pipeline_analyze:doc_analyze_streaming:183 - Pipeline processing-window multi-file run. doc_count=1, total_pages=12, window_size=64, total_batches=1
2026-05-09 09:12:58.938 | INFO     | mineru.backend.pipeline.pipeline_analyze:doc_analyze_streaming:235 - Pipeline processing window batch 1/1: 12/12 pages, batch_pages=12, doc_slices=doc0:1-12
2026-05-09 09:12:58.939 | INFO     | mineru.backend.pipeline.pipeline_analyze:batch_image_analyze:328 - GPU Memory: 1 GB, Batch Ratio: 1. 
2026-05-09 09:12:58.939 | INFO     | mineru.backend.pipeline.model_init:__init__:207 - DocAnalysis init, this may take some times......
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:12:59,851 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:00,845 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:01,841 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:02,767 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:04,033 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:04,919 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:05,979 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:06,825 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:07,760 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:08,620 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:09,448 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:13:10,450 - modelscope - INFO - Target directory already exists, skipping creation.
2026-05-09 09:13:10.470 | INFO     | mineru.backend.pipeline.model_init:__init__:260 - DocAnalysis init done!
2026-05-09 09:13:10.470 | INFO     | mineru.backend.pipeline.pipeline_analyze:custom_model_init:83 - model init cost: 11.531244039535522
Layout Predict: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 12/12 [00:02<00:00,  4.16it/s]
MFR Predict: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████| 102/102 [00:47<00:00,  2.15it/s]
Table-ocr det: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [00:00<00:00, 24.41it/s]
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:14:02,209 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:14:03,144 - modelscope - INFO - Target directory already exists, skipping creation.
Table-ocr rec ch: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████| 99/99 [00:00<00:00, 99.01it/s]
Table-wireless Predict: 100%|███████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [00:00<00:00, 14.78it/s]
Table-wired Predict:   0%|                                                                                                          | 0/1 [00:00<?, ?it/s]Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:14:05,469 - modelscope - INFO - Target directory already exists, skipping creation.
Table-wired Predict: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████| 1/1 [00:01<00:00,  1.21s/it]
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:14:06,604 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-09 09:14:07,692 - modelscope - INFO - Target directory already exists, skipping creation.
OCR-det ch: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 42/42 [00:03<00:00, 12.36it/s]
Seal Predict: 0it [00:00, ?it/s]
OCR-rec Predict: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████| 57/57 [00:00<00:00, 93.73it/s]
Processing pages: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████| 12/12 [00:01<00:00,  7.93it/s]
2026-05-09 09:14:14.019 | INFO     | mineru.cli.client:run_planned_task:807 - Completed batch 1/1 | Processed 12/12 pages | 1 of 1 batch finished | task#1 [Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En]
INFO:     Shutting down
INFO:     Waiting for application shutdown.
INFO:     Application shutdown complete.
INFO:     Finished server process [2957300]
Done.

```

the output is not merely a file, but a directory containing the parsed markdown file and the original PDF file, which is useful for you to check the parsing quality by comparing the markdown file with the original PDF file.

Its hierarchy is as follows:

```bash
# test history
./test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En
└── auto
    ├── images
    │   ├── 108ab5199f55198dabe5235a25c47d5948d7a1f94c7f8ad21820772ea5f302e4.jpg
    |   ├── # lots of .jpg files
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_content_list.json
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_content_list_v2.json
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_layout.pdf
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.md
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_middle.json
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_model.json
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_origin.pdf
    └── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_span.pdf

3 directories, 45 files

```

And we only use the `.md files and _content_list_v2.json/_content_list.json files` for further processing like structuring.

So, if you only do not want to spare time dealing with the rest of the files, you can use the `--clear` argument to strip anything unnecessary.

you can run the following command 
```bash
paperflow pdf-parse -i ./test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.pdf  -o ./test/Other_database/  --clear
```

the log shows the same
```bash
Running: mineru -p /data2/pyPaperFlow/test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.pdf -o /data2/pyPaperFlow/test/Other_database -b pipeline
2026-05-11 14:19:36.551 | INFO     | mineru.cli.client:run_orchestrated_cli:874 - Started local mineru-api at http://127.0.0.1:51213
2026-05-11 14:19:37.548 | INFO     | __main__:create_app:260 - Request concurrency limited to 3
Start MinerU FastAPI Service: http://127.0.0.1:51213
API documentation: http://127.0.0.1:51213/docs
INFO:     Started server process [3574652]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://127.0.0.1:51213 (Press CTRL+C to quit)
2026-05-11 14:19:38.557 | INFO     | mineru.cli.client:run_planned_task:771 - Submitting batch 1/1 | 1 document, 12 pages in this batch | 12 pages total | task#1 [Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En]
2026-05-11 14:19:39.860 | INFO     | mineru.backend.pipeline.pipeline_analyze:doc_analyze_streaming:183 - Pipeline processing-window multi-file run. doc_count=1, total_pages=12, window_size=64, total_batches=1
2026-05-11 14:19:41.235 | INFO     | mineru.backend.pipeline.pipeline_analyze:doc_analyze_streaming:235 - Pipeline processing window batch 1/1: 12/12 pages, batch_pages=12, doc_slices=doc0:1-12
2026-05-11 14:19:41.236 | INFO     | mineru.backend.pipeline.pipeline_analyze:batch_image_analyze:328 - GPU Memory: 1 GB, Batch Ratio: 1. 
2026-05-11 14:19:41.236 | INFO     | mineru.backend.pipeline.model_init:__init__:207 - DocAnalysis init, this may take some times......
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:42,150 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:43,195 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:44,141 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:44,949 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:46,272 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:47,155 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:48,168 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:48,942 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:49,836 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:50,785 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:51,640 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:19:52,685 - modelscope - INFO - Target directory already exists, skipping creation.
2026-05-11 14:19:52.700 | INFO     | mineru.backend.pipeline.model_init:__init__:260 - DocAnalysis init done!
2026-05-11 14:19:52.700 | INFO     | mineru.backend.pipeline.pipeline_analyze:custom_model_init:83 - model init cost: 11.464261293411255
Layout Predict: 100%|██████████████████████████████████| 12/12 [00:03<00:00,  3.96it/s]
MFR Predict: 100%|███████████████████████████████████| 102/102 [00:48<00:00,  2.10it/s]
Table-ocr det: 100%|█████████████████████████████████████| 2/2 [00:00<00:00, 24.59it/s]
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:20:45,866 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:20:46,730 - modelscope - INFO - Target directory already exists, skipping creation.
Table-ocr rec ch: 100%|███████████████████████████████| 99/99 [00:00<00:00, 102.41it/s]
Table-wireless Predict: 100%|████████████████████████████| 2/2 [00:00<00:00, 19.19it/s]
Table-wired Predict:   0%|                                       | 0/1 [00:00<?, ?it/s]Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:20:48,937 - modelscope - INFO - Target directory already exists, skipping creation.
Table-wired Predict: 100%|███████████████████████████████| 1/1 [00:01<00:00,  1.11s/it]
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:20:49,964 - modelscope - INFO - Target directory already exists, skipping creation.
Downloading Model from https://www.modelscope.cn to directory: /home/nicai_zht/.cache/modelscope/hub/models/OpenDataLab/PDF-Extract-Kit-1.0
2026-05-11 14:20:50,892 - modelscope - INFO - Target directory already exists, skipping creation.
OCR-det ch: 100%|██████████████████████████████████████| 42/42 [00:03<00:00, 12.27it/s]
Seal Predict: 0it [00:00, ?it/s]
OCR-rec Predict: 100%|█████████████████████████████████| 57/57 [00:00<00:00, 98.48it/s]
Processing pages: 100%|████████████████████████████████| 12/12 [00:01<00:00,  8.16it/s]
2026-05-11 14:20:57.180 | INFO     | mineru.cli.client:run_planned_task:807 - Completed batch 1/1 | Processed 12/12 pages | 1 of 1 batch finished | task#1 [Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En]
INFO:     Shutting down
INFO:     Waiting for application shutdown.
INFO:     Application shutdown complete.
INFO:     Finished server process [3574652]
Done.
✅Removed 5 source files. Only .md and necessary .json files are kept in the output directory test/Other_database.

```

now the output directory is much cleaner, it will be much more disk-space-saving for you to perform batch processing on a large number of papers.
```bash
Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En
└── auto
    ├── images
    │   ├── 108ab5199f55198dabe5235a25c47d5948d7a1f94c7f8ad21820772ea5f302e4.jpg
    │   ├── # lots of .jpg files
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_content_list.json
    ├── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_content_list_v2.json
    └── Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En.md

3 directories, 40 files
```

Because the JSON file is hard to parse into a Section-based hierarchical structure markdown file, all we need is the markdown file. But you can use other original output files for analysis or debugging.

Remember: 
- We only need the structured markdown file like what we have done in PMC papers.
Structured section like `Introduction`, `Methods`, `Results`, `Discussion`, `Conclusion`. 
- All begin with a JSON file, we parse everything and do post-processing job only for a JSON output, which contains the metadata 
and the content sections we mentioned above.
- We do selection/aggregation step ONLY on the final JSON file mentioned above, and you provide a YAML configuration file to specify which sections to extract and compile into the final markdown file. 


Because the `original JSON file generated by mineru` is hard to use directly, so we currently only use the markdown and _content_list_v2.json files for further structuring.

## 🧬 Case 7: Parse MinerU content_list_v2.json into canonical sectioned JSON

### Overview / 概述

After running MinerU's `pdf-parse` (Case 6), you get a `content_list_v2.json` file. This JSON contains raw, page-by-page block data from the PDF — titles, paragraphs, images, tables, etc. — but with no semantic structure. Section headings like "1. Introduction" or "Experimental Section" are just text strings.

The `mineru-parse` command transforms this flat JSON into a structured, canonical JSON where every section is classified into a standard academic type (abstract, introduction, methods, results, discussion, etc.). Metadata (title, authors, year, DOI, journal) and figure captions are also extracted.

   
---

    


A sliding cursor tracks document order to reduce false matches (a late "Methods" heading in an unusual position is still recognized).



---
   
### Usage / 使用方法      

```bash
# Regex backend (default, no API key needed)
paperflow mineru-parse -i content_list_v2.json -o paper.json

# AI backend with Claude
export ANTHROPIC_API_KEY="sk-ant-..."
paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai

# AI backend with DeepSeek    
export OPENAI_API_KEY="sk-..."
paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai \
  --base-url https://api.deepseek.com --model deepseek-v4-pro

# AI backend with university proxy / 大学代理
paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai \
  --base-url https://models.sjtu.edu.cn/api/v1 --model deepseek-chat --api-key your-key

# Custom config / 使用自定义配置
paperflow mineru-parse -i content_list_v2.json -o paper.json --config my_rules.yaml
```


**15 canonical section types / 15 种规范章节类型：**

`abstract` `introduction` `results` `discussion` `methods` `conclusion` `supplementary` `availability` `funding` `acknowledgements` `author_contributions` `keywords` `conflicts` `references` `other`

---  
 

### Notes & Limitations / 注意事项与局限


目前建议是使用自己测试的样本不断完善结构解析的边界情况，直到你觉得大部分的论文都能被正确解析了，再进行批量处理。



the following table shows


| CLI | Log |
| :--- | :--- |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En/auto/Zhu_2025_AdvancedScience_Accurate_Generation_of_Conformational_En_content_list_v2.json -o idpfold1.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml` | Using regex backend with configurable aliases<br>Parsed 10 sections -> idpfold1.json<br>Sections: abstract(Abstract), introduction(Introduction), results(Results), discussion(Discussion), methods(Methods), supplementary(Supplementary Material), availability(Data & Code Availability), acknowledgements(Acknowledgements), keywords(Keywords), conflicts(Competing Interests) |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/staring/auto/staring_content_list_v2.json  -o staring.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml`  | Using regex backend with configurable aliases<br>Parsed 14 sections -> staring.json<br>Sections: abstract(Abstract), discussion(Discussion), methods(Methods), supplementary(Supplementary Material), availability(Data & Code Availability), other(Online content), other(Article), other(Additional information), other(Article), other(Statistics), other(Software and code), other(Data), other(Field-specific reporting), other(Plants)  |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/idpfold2/auto/idpfold2_content_list_v2.json  -o idpfold2.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml` | Using regex backend with configurable aliases<br>Parsed 14 sections -> idpfold2.json<br>Sections: abstract(Abstract), introduction(Introduction), discussion(Discussion), methods(Methods), availability(Data & Code Availability), acknowledgements(Acknowledgements), keywords(Keywords), conflicts(Competing Interests), references(References), other(Overview), other(Predicting global compaction across the order-disorder continuum), other(Fitting global and local experimental observations), other(Modelling multiple conformations for protein assemblies), other(Prediction conformational changes in IDR-binding) |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/disobind/auto/disobind_content_list_v2.json  -o disobind.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml`  | Using regex backend with configurable aliases<br>Parsed 36 sections -> disobind.json<br>Sections: abstract(Abstract), introduction(Introduction), discussion(Discussion), methods(Methods), supplementary(Supplementary Material), availability(Data & Code Availability), funding(Funding), acknowledgements(Acknowledgements), author_contributions(Author Contributions), conflicts(Competing Interests), references(References), other(Inter-protein contact map prediction), other(Interface residue prediction), other(Coarse-graining improves the performance), other(Comparison to AlphaFold2 and AlphaFold3), other(Using different ipTM cutoffs for AF2 and AF3), other(AF2 performs better than AF3), other(Combining Disobind and AlphaFold2 predictions), other(Performance by residue types), other(Comparison with interface predictors for IDRs), other(Protein language models allow for a shallow architecture), other(Diversity and inclusion statement), other(Declaration of generative AI and AI-assisted technologies), other(Tables), other(Key resources table), other(Gathering PDB structures of IDRs in complexes), other(Defining binary complexes containing IDRs), other(Creating merged binary complexes), other(Notations), other(Inputs and outputs for training), other(Projection dimension), other(Number of layers in the MLP), other(,  for SE loss), other(Disobind+AF2 predictions), other(Performance by residue type), other(Comparison with interface predictors for IDRs)  |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/alphafold3/auto/alphafold3_content_list_v2.json  -o alphafold3.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml`   |  Using regex backend with configurable aliases<br>Parsed 14 sections -> alphafold3.json<br>Sections: abstract(Abstract), discussion(Discussion), methods(Methods), supplementary(Supplementary Material), availability(Data & Code Availability), other(Model limitations), other(Online content), other(Metrics), other(Nucleic acid prediction baseline), other(Model performance analysis and visualization), other(Additional information), other(Article), other(Statistics), other(Software and code)  |
| `paperflow mineru-parse -i /data2/pyPaperFlow/test/Other_database/2409.02240v1/auto/2409.02240v1_content_list_v2.json -o 2409.02240v1.json --backend regex --config /data2/pyPaperFlow/src/pyPaperFlow/integrations/mineru_config.yaml` | Using regex backend with configurable aliases<br>Parsed 13 sections -> 2409.02240v1.json<br>Sections: abstract(Abstract), introduction(Introduction), methods(Methods), conclusion(Conclusion), acknowledgements(Acknowledgements), author_contributions(Author Contributions), conflicts(Competing Interests), references(References), other(Ensemble Generation), other(Ensemble Validation), other(Experimental Observables), other(ML Ensemble Generation and Validation), other(Software and Data Repositories for IDPs/IDRs)  |

you can check all the output JSON files here : [parse](../test/Other_database/parse/)



### 4. 搜索并获取 arXiv 论文
如果你只想先拿到 ID，可以先搜索；如果想同时获取元数据和 PDF，可以直接 fetch。

```bash
paperflow arxiv-search "deep learning for biology" --max-results 10
paperflow arxiv-fetch "deep learning for biology" --max-results 10 --download-pdf
paperflow arxiv-fetch "deep learning for biology" --max-results 10 --download-pdf --backend paperscraper
```

常用参数：

- `--start-date` / `--end-date`：按 `YYYY-MM-DD` 格式限制日期范围。
- `--backend`：可选 `native`（内置的 httpx 方案）或 `paperscraper`（安装了第三方包时可用）。
- `--output-dir`：把 ID 列表或抓取结果保存到其他目录。
- `--no-download-pdf`：只保存元数据，不下载 PDF。

日期过滤示例：

```bash
paperflow arxiv-fetch "protein folding" --start-date 2024-01-01 --end-date 2024-12-31 -o ./papers/arxiv
```

搜索结果会保存为 `searched_arxiv_ids.txt`。抓取结果会按 `source/year/source_id/` 结构保存，包含 JSON 元数据，PDF 则按可用情况尽量下载。

# 🌟🌟🌟
#### arXiv 命令变体与使用示例

- **arxiv-search**: 仅检索匹配的 arXiv 记录并输出 ID 列表（不下载内容）。

    用法示例：

    ```bash
    paperflow arxiv-search "protein folding" --max-results 50 --start-date 2024-01-01 --end-date 2024-12-31
    # 将会在默认存储目录下生成 searched_arxiv_ids.txt，或使用 --output-dir 指定保存位置
    ```
    
- **arxiv-fetch**: 检索并保存每篇论文的标准化元数据（JSON），可选地下载 PDF 文件（默认开启）。

    常用选项：
    - `--download-pdf/--no-download-pdf`：是否下载 PDF（默认 `--download-pdf`）。
    - `--backend`：`native`（默认，使用 arXiv Atom API）或 `paperscraper`（需安装 `paperscraper` 包）。
    - `--output-dir`：指定保存结果的目录（默认使用全局存储目录）。
    - `--start-date` / `--end-date`：按 `YYYY-MM-DD` 限制提交时间范围。

    用法示例：

    ```bash
    # 仅保存元数据（不下载 PDF）
    paperflow arxiv-fetch "deep learning for biology" --max-results 20 --no-download-pdf -o ./papers/arxiv

    # 使用 paperscraper 后端并下载 PDF
    paperflow arxiv-fetch "deep learning for biology" --max-results 20 --download-pdf --backend paperscraper -o ./papers/arxiv
    ```

- **输出与存储**：
    - 元数据：每篇论文保存为 `{source_id}.json`，包含 `title`, `authors`, `abstract`, `published_date`, `landing_url`, `pdf_url` 等字段（存储路径示例：`{output_dir}/arxiv/2024/2301.01234v1/2301.01234v1.json`）。
    - PDF：如果可用且下载成功，则保存为 `{source_id}.pdf`，并在对应 JSON 中更新 `pdf_downloaded` 和 `pdf_path` 字段。

- **注意事项**：
    - arXiv 的抓取流程只负责元数据标准化与 PDF 下载；当前仓库没有内建将 arXiv PDF 自动解析为 Markdown/结构化全文的步骤。若需后续文本解析，请在下载后接入 PDF 解析器（例如 `pdfplumber`、`minerU`、或 OCR/布局解析管线），并将解析结果保存为 `*_parsed.md` 或结构化 JSON，以便 `merge` 等下游工具使用。


### 5. 搜索并获取 bioRxiv 论文
bioRxiv 目前走 Crossref（openRxiv）服务端直接检索，不再先拉取大范围日期窗口再在本地做匹配。

```bash
paperflow biorxiv-search "AlphaFold AND structure" --max-results 10
paperflow biorxiv-fetch "AlphaFold AND structure" --start-date 2026-01-01 --end-date 2026-01-31 --download-pdf
```

常用参数：

- `--start-date` / `--end-date`：按 `YYYY-MM-DD` 格式限制日期范围。
- `--output-dir`：把 ID 列表或抓取结果保存到其他目录。
- `--no-download-pdf`：只保存元数据，不下载 PDF。

兼容性说明：

- `--window-days` 作为 CLI 兼容参数保留，但当前 Crossref 检索路径不会使用该参数。

示例：  

```bash
paperflow biorxiv-fetch "protein interaction" --max-results 50 -o ./papers/biorxiv
```

搜索结果会保存为 `searched_biorxiv_ids.txt`。抓取结果会按 `source/year/source_id/` 结构保存，包含 JSON 元数据，并在可用时下载 PDF。




### 4. Search and Fetch arXiv Papers
Search arXiv first if you only want IDs, or fetch metadata and PDFs in one step.

```bash
paperflow arxiv-search "deep learning for biology" --max-results 10
paperflow arxiv-fetch "deep learning for biology" --max-results 10 --download-pdf
paperflow arxiv-fetch "deep learning for biology" --max-results 10 --download-pdf --backend paperscraper
```

Useful options:

- `--start-date` and `--end-date`: limit results to a date window in `YYYY-MM-DD` format.
- `--backend`: choose `native` for the built-in httpx-backed arXiv API path, or `paperscraper` to use the optional third-party adapter when installed.
- `--output-dir`: save the ID list or fetched records to a different directory.
- `--no-download-pdf`: skip PDF download and save metadata only.

Example with a date filter:

```bash
paperflow arxiv-fetch "protein folding" --start-date 2024-01-01 --end-date 2024-12-31 -o ./papers/arxiv
```

Search output is saved as `searched_arxiv_ids.txt`. Fetched records are stored under `source/year/source_id/` with JSON metadata and, when available, a PDF copy.

### 5. Search and Fetch bioRxiv Papers
bioRxiv now uses direct server-side query via Crossref (openRxiv records), rather than pulling large date windows first and filtering locally.

```bash
paperflow biorxiv-search "AlphaFold AND structure" --max-results 10
paperflow biorxiv-fetch "AlphaFold AND structure" --start-date 2026-01-01 --end-date 2026-01-31 --download-pdf
```

Useful options:

- `--start-date` and `--end-date`: limit results to a date window in `YYYY-MM-DD` format.
- `--output-dir`: save the ID list or fetched records to a different directory.
- `--no-download-pdf`: skip PDF download and save metadata only.

Compatibility note:

- `--window-days` is kept for CLI compatibility but is not used by the current Crossref-backed bioRxiv search path.

Example:

```bash
paperflow biorxiv-fetch "protein interaction" --max-results 50 -o ./papers/biorxiv
```

Search output is saved as `searched_biorxiv_ids.txt`. Fetched records are stored under `source/year/source_id/` with JSON metadata and, when available, a PDF copy.