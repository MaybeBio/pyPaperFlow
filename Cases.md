# Test Cases

[Back to README](./README.md)

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

![alt text](./figs/pubmed.png)

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

you can check the result here: [alphafold3_ensemble_meta](./test/alphafold3_ensemble_meta/)

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

Hi @MaybeBio 

You should activate proxy when you using OpenAI API format.
> back to main menu and press "P"

If It still not work, you can type the CLI "cc-switch proxy show", and find the field
`
Issue #49 manual Claude setup:
- ANTHROPIC_BASE_URL=http://127.0.0.1:15721
- ANTHROPIC_AUTH_TOKEN=proxy-placeholder
- Keep the real upstream base URL and API key in the selected Claude provider inside cc-switch.
`

set ANTHROPIC_BASE_URL and ANTHROPIC_AUTH_TOKEN according to the output.

By the way, I recommend you use OpenAI response API instead of completion API. Because the function call has some problem in completion API. I've already fix this issue, and it should be involved in the next version.