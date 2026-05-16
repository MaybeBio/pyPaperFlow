---
name: pubmed-query-builder
description: Construct optimized PubMed search queries through systematic iterative refinement. Use when user mentions literature search, building PubMed queries, finding papers, constructing search strings, or needs systematic literature review methodology. Trigger for requests like "search for X papers", "build a query for Y", "find literature on Z", or any biomedical literature search needs.
---

# PubMed Query Builder

Systematic construction of PubMed search queries through iterative refinement and validation.

## Iron Laws

1. **Requirement Elicitation Before Construction**: Never accept vague requests like "find cancer papers". Must systematically elicit specific research dimensions: disease, intervention, outcomes, population, study design.

2. **MeSH Vocabulary Verification Mandatory**: For ALL biomedical and clinical concepts, user must verify terms against the MeSH Browser (https://meshb.nlm.nih.gov/). NEVER proceed with unverified keywords. Use Entry Terms (synonyms) identified in MeSH records to expand search coverage.

3. **Multi-Dimensional Query Construction**: Comprehensive queries must integrate: verified MeSH terms (including Entry Terms), keywords for non-MeSH concepts, conceptual hierarchies (using [noexp] appropriately), temporal constraints, publication types, population restrictions, and exclusion criteria.

4. **Syntax Validation**: Final output must conform to PubMed syntax specifications: correct field tags, proper Boolean operator nesting, valid date formats, appropriate use of quotation marks and wildcards.

## Hard Gates

Before proceeding to query construction, one of the following must be satisfied:

- ✓ User provides a specific research question or hypothesis
- ✓ User identifies a specific disease/phenotype name
- ✓ User specifies a molecular entity (gene, protein, drug)
- ✓ User names a specific research methodology or technique

**MANDATORY MeSH Verification**:
Before ANY query construction, user must have verified their concepts at https://meshb.nlm.nih.gov/. If user has not done this, construction CANNOT proceed. This is non-negotiable.

## Rationalization Table

| Potential Thought | Actual Reality |
|-------------------|----------------|
| "User only gave a keyword, simple search suffices" | Inadequate queries yield irrelevant results, wasting user time and potentially missing critical literature |
| "MeSH terms are complex, keywords work fine" | Keyword-only searches miss relevant articles and retrieve irrelevant ones, particularly with common terms or new concepts |
| "I can guess the MeSH term from the keyword" | Incorrect MeSH terms yield fundamentally flawed searches; verification at https://meshb.nlm.nih.gov/ is mandatory |
| "Entry Terms are unnecessary, just use the MeSH Heading" | Entry Terms capture synonyms and alternative terminology used by authors; omitting them reduces recall |
| "Complex syntax confuses users, provide simple version" | Users need precise results, not simplified syntax; accuracy outweighs perceived complexity |

## Red Flags

Immediately halt when these thoughts arise:

- "This concept is too specific, no MeSH exists" → Search broader MeSH terms and apply restrictions; verify at https://meshb.nlm.nih.gov/
- "User hasn't verified MeSH terms, proceed anyway" → STOP. Verification is mandatory before any construction
- "Entry Terms are optional, skip them" → Entry Terms are essential for comprehensive retrieval; include ALL listed Entry Terms
- "This query is sufficiently complex, no further optimization needed" → Complexity is not the metric; precision and recall balance determines quality
- "Users cannot understand advanced syntax, provide basic version" → Users require accurate results, provide both syntax and explanation

## PubMed Search Syntax Reference

### Field Tags

| Tag | Field | Description |
|-----|-------|-------------|
| `[mh]` | MeSH Terms | NLM Medical Subject Headings controlled vocabulary |
| `[majr]` | MeSH Major Topic | Articles where term is major focus |
| `[sh]` | MeSH Subheadings | Qualifiers attached to MeSH terms |
| `[ti]` | Title | Article title words |
| `[ab]` | Abstract | Abstract words |
| `[tiab]` | Title/Abstract | Combined title and abstract search |
| `[tw]` | Text Words | All fields including MeSH, Publication Types, Substance Names |
| `[au]` | Author | Author names (format: lastname initials) |
| `[ta]` | Journal | Journal title or abbreviation |
| `[dp]` | Publication Date | Date format: YYYY/MM/DD |
| `[pt]` | Publication Type | Article type (e.g., Clinical Trial, Review) |
| `[la]` | Language | Publication language |
| `[sb]` | Subset | Subject, citation status, or journal category filters |

### Boolean Operators

- **AND**: Retrieves results containing ALL search terms (default between concepts)
- **OR**: Retrieves results containing AT LEAST ONE search term
- **NOT**: Excludes terms from search results

PubMed processes searches left-to-right; use parentheses for nesting.

### MeSH Terms and Hierarchies

MeSH (Medical Subject Headings) is NLM's controlled vocabulary thesaurus for indexing PubMed citations. MeSH terms are arranged hierarchically with specific terms beneath broader terms.

**Automatic Expansion**: By default, MeSH searches include all more specific terms in the hierarchy.

**Disabling Expansion**: Use `[noexp]` suffix to exclude automatic inclusion of indented terms:
```
neoplasms[mh:noexp]  # Searches only "neoplasms", not specific cancer types
```

**Exact Match Only**: For exact term matching without mapping to related terms, use double quotes:
```
"neoplasms"[mh]
```

**MeSH with Subheadings**: Direct attachment using forward slash:
```
hypertension/therapy[mh]  # Hypertension with therapy subheading
```

Subheading abbreviations may be used:
```
hypertension/th[mh]  # "th" = therapy
```

**Major MeSH Topic**: Use `[majr]` to retrieve articles where the MeSH term is a major focus:
```
hypertension[majr]
```

### Phrase Searching

**Exact Phrase**: Enclose in double quotes:
```
"kidney allograft"
```

Note: If phrase is not in the phrase index, PubMed ignores quotes and processes using Automatic Term Mapping.

**Phrase with Field Tag**:
```
"kidney allograft"[tiab]
```

**Bypassing Phrase Index**: Use proximity search with distance 0 for phrases not in phrase index:
```
"cognitive impairment in multiple sclerosis"[tiab:~0]
```

### Wildcards and Truncation

Use asterisk (*) to substitute for 0 or more characters:

```
vaccin*  # Matches: vaccine, vaccination, vaccinating, etc.
organiz*ion*  # Matches: organization, organisation, etc.
```

**Requirements**:
- Terms must begin with at least 3-4 characters before wildcard
- Multiple wildcards allowed in single term

**Phrase with Wildcards**:
```
"breast feed*"
"colo* cancer*"
```

Wildcards disable Automatic Term Mapping and MeSH hierarchical expansion.

### Date Searching

**Single Date**:
```
1998/03/06[dp]
1998[dp]  # Month and day optional
1998/03[dp]
```

**Date Range**:
```
1996:1998[dp]
1998/01:1998/04[dp]
```

**Relative Dates**:
```
"last 5 years"[dp]
"last 30 days"[dp]
"last 12 months"[dp]
```

Date ranges include both print and electronic publication dates.

### Publication Types

Common publication types include:

| Publication Type | Tag |
|-----------------|-----|
| Clinical Trial | `clinical trial[pt]` |
| Randomized Controlled Trial | `randomized controlled trial[pt]` |
| Systematic Review | `systematic review[pt]` |
| Meta-Analysis | `meta-analysis[pt]` |
| Review | `review[pt]` |
| Case Report | `case reports[pt]` |
| Guideline | `practice guideline[pt]` |

Publication types are hierarchical; use `[noexp]` to disable automatic inclusion of specific types:
```
clinical trial[pt:noexp]
```

### Character Conversions

Special characters in searches:
- `()` - Boolean nesting
- `[]` - Field tag qualification
- `&` - Boolean AND
- `|` - Boolean OR
- `/` - MeSH/Subheading combinations
- `:` - Range operation
- `""` - Phrase search
- `#` - History statement (when followed by number)
- `*` - Wildcard

Characters converted to spaces: `!`, `#`, `$`, `%`, `+`, `-`, `.`, `,`, `;`, `<`, `>`, `=`, `?`, `\`, `^`, `_`, `{}`, `~`, `'`

### Filters and Subsets

**Systematic Review Filter**:
```
(systematic[sb] OR systematic review[pt])
```

**Full Text Filter**:
```
"full text"[sb]
```

**Free Full Text Filter**:
```
"free full text"[sb]
```

**MEDLINE Subset**:
```
medline[sb]
```

**Exclude Preprints**:
```
NOT preprint[pt]
```

### Clinical Study Categories

**Therapy (Sensitive/Broad)**:
```
((clinical[ti] AND trial[ti]) OR clinical trials as topic[mh] OR clinical trial[pt] OR random*[ti] OR random allocation[mh] OR therapeutic use[sh])
```

**Therapy (Specific/Narrow)**:
```
(randomized controlled trial[pt] OR (randomized[ti] AND controlled[ti] AND trial[ti]))
```

## Systematic Query Construction Workflow

### Phase 1: Requirement Elicitation

Employ systematic questioning to clarify user needs:

**Research Question**
- "What is your specific research question or hypothesis?"
- "What knowledge gap are you attempting to address?"

**Core Concept Identification**
- "What are the primary biomedical concepts involved?"
- "Identify diseases/conditions, interventions/exposures, outcomes, population characteristics."

**CRITICAL: MeSH Term Verification Instructions**
- For each concept identified, you MUST verify the official MeSH terminology
- Go to https://meshb.nlm.nih.gov/ and search for each concept
- Do NOT assume your keyword is the correct MeSH term
- Record the MeSH Heading, Entry Terms (all synonyms), and Unique ID
- Use this verified information to construct the search terms

**Example Verification Process**:
- User concept: "disordered protein"
- MeSH Browser search: https://meshb.nlm.nih.gov/ → search "disordered protein"
- Found: "Intrinsically Disordered Proteins" (D064267)
- Entry Terms: "Intrinsically Disordered Protein", "Natively Unfolded Protein", "Natively Unfolded Proteins", "Unstructured Protein", "Unstructured Proteins"
- Query construction: Use official MeSH Heading AND all Entry Terms

**Dimensional Specification**
- Population: "What are the demographic and clinical characteristics of interest?" (age, sex, disease stage, comorbidities)
- Intervention/Exposure: "What specific interventions, drugs, procedures, or exposures?"
- Comparison: "Is there a comparator or control condition?"
- Outcomes: "What endpoints or outcomes are of interest?"
- Study Design: "What study methodologies are most relevant?" (RCT, cohort, case-control, systematic review, meta-analysis)
- Temporal Scope: "What publication date range is appropriate?"
- Language: "Are language restrictions required?"
- Exclusions: "What types of studies should be excluded?"

### Phase 2: Concept Mapping and Term Identification

For each core concept:

**1. MeSH Vocabulary Verification (MANDATORY)**
- User MUST verify each concept at https://meshb.nlm.nih.gov/
- Search by concept name to find the official MeSH Heading
- Record the MeSH Heading (preferred term), Unique ID (UI), Tree Number(s)
- Document ALL Entry Terms listed in the MeSH record
- Note Scope Note for precise definition
- Verify Year Introduced to assess coverage period
- Use Entry Terms as OR-connected synonyms in the query

**Example: "disordered protein"**
- MeSH Browser search reveals: "Intrinsically Disordered Proteins" (D064267)
- Entry Terms include: "Intrinsically Disordered Protein", "Natively Unfolded Protein", "Natively Unfolded Proteins", "Unstructured Protein", "Unstructured Proteins"
- Query component: `"Intrinsically Disordered Proteins"[mh] OR "Intrinsically Disordered Protein"[mh] OR "Natively Unfolded Protein"[mh] OR "Natively Unfolded Proteins"[mh] OR "Unstructured Protein"[mh] OR "Unstructured Proteins"[mh]`

**2. Keyword Strategy for Non-MeSH Concepts**
- Identify terms without MeSH equivalents (new concepts, brand names, proprietary acronyms)
- Include variants, abbreviations, and spelling variations
- Consider recent terminology not yet incorporated into MeSH
- Use field tags [tiab] or [tw] for these terms

**3. Field Tag Selection**
- `[mh]` for verified MeSH terms (maximum precision)
- `[majr]` for major MeSH topics (articles where concept is primary focus)
- `[tiab]` for keywords in title/abstract (balanced precision/recall)
- `[tw]` for comprehensive text word searching (maximum recall, lower precision)
- `[pt]` for publication types (study design methodology)

### Phase 3: Query Structure Construction

Employ standardized PICO-based structure:

```
(Population OR Demographics) AND (Intervention OR Exposure) AND (Comparator) AND (Outcome) AND (Study Design Filters)
```

**Syntax Principles**:
1. Enclose each concept group in parentheses
2. Use uppercase Boolean operators (AND, OR, NOT)
3. OR connects synonyms within a concept group
4. AND connects distinct concept groups
5. Place NOT operators at the end of concept groups
6. Apply temporal and design filters after core concepts

**Example Structure**:
```
(("hypertension"[mh] OR "high blood pressure"[tiab]) AND 
("ACE inhibitors"[mh] OR "angiotensin converting enzyme inhibitors"[tiab] OR 
("lisinopril"[tw] OR "enalapril"[tw] OR "ramipril"[tw])) AND 
("cardiovascular mortality"[tiab] OR "myocardial infarction"[tiab] OR "stroke"[tiab]) AND 
("randomized controlled trial"[pt] OR "clinical trial"[pt]) AND 
2015/01/01:2024/12/31[dp])
```

### Phase 4: Iterative Refinement

Present initial query structure and elicit feedback:

**Relevance Assessment**
- "Does this query accurately represent your research intent?"
- "Are there concepts that should be added or removed?"

**Scope Optimization**
- "Should the search be expanded (more terms, broader MeSH) or narrowed (more specific filters)?"
- "Are temporal constraints appropriate?"
- "Should publication type filters be added or modified?"

**Exclusion Criteria**
- "What types of studies or populations should be excluded?"

**Validation Approach**
- "Will you be validating the query against a set of known relevant articles?"
- "Do you have a gold standard article set for testing?"

Iteratively modify based on user feedback until consensus is achieved.

## Output Format

Provide comprehensive documentation:

```markdown
## PubMed Query

```
(final_query_string)
```

## Query Decomposition

| Concept Component | Search Terms | Field Tags | Rationale |
|-------------------|--------------|------------|-----------|
| Disease/Condition | "hypertension"[mh] OR "high blood pressure"[tiab] | MeSH, TIAB | MeSH ensures standardized terminology; TIAB captures variants |
| Intervention | "ACE inhibitors"[mh] OR ("lisinopril"[tw] OR "enalapril"[tw]) | MeSH, TW | MeSH class term plus individual drug names |
| Outcome | "cardiovascular mortality"[tiab] OR "stroke"[tiab] | TIAB | Clinical endpoints expressed in author terminology |
| Study Design | "randomized controlled trial"[pt] OR "clinical trial"[pt] | PT | Restricts to clinical trial methodologies |
| Temporal Range | 2015/01/01:2024/12/31[dp] | DP | Limits to recent decade |

## PubMed Search URL

https://pubmed.ncbi.nlm.nih.gov/?term=(URL_encoded_query)

## Anticipated Results

Estimated range: X-Y citations based on concept specificity and applied filters.

## Optimization Recommendations

[Optional suggestions for further refinement based on initial results or known literature]
```

## Validation Checklist

Before declaring completion, verify:

- [ ] User has verified ALL core concepts at https://meshb.nlm.nih.gov/
- [ ] Each concept uses the official MeSH Heading
- [ ] ALL Entry Terms listed in MeSH records are included
- [ ] MeSH Unique IDs are documented for reference
- [ ] All field tags are correctly applied ([mh], [majr], [tiab], [tw], [pt], [dp], etc.)
- [ ] Boolean operators are uppercase and properly nested
- [ ] Parentheses correctly encapsulate concept groups
- [ ] Exclusion terms (NOT) are positioned appropriately
- [ ] Date ranges follow YYYY/MM/DD:YYYY/MM/DD format
- [ ] Publication types conform to PubMed conventions
- [ ] User confirms query accurately represents research intent
- [ ] Search URL is functional and properly encoded
- [ ] Documentation provides rationale for each decision

## Reference Resources

**MeSH Database**: https://www.ncbi.nlm.nih.gov/mesh/

**PubMed Advanced Search Builder**: https://pubmed.ncbi.nlm.nih.gov/advanced/

**PubMed Clinical Queries**: https://pubmed.ncbi.nlm.nih.gov/clinical/

**Single Citation Matcher**: https://pubmed.ncbi.nlm.nih.gov/citmatch/

**PubMed Help Documentation**: https://pubmed.ncbi.nlm.nih.gov/help/

**E-utilities API**: https://www.ncbi.nlm.nih.gov/books/n/helpeutils/chapter1/

**MeSH Browser**: https://meshb.nlm.nih.gov/

**Publication Types**: https://www.nlm.nih.gov/mesh/pubtypes.html

**Clinical Study Categories**: https://www.ncbi.nlm.nih.gov/books/n/helpeutils/chapter2/#Clinical_Queries_Filters