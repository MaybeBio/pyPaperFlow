**来源页面:** [Help - PubMed](https://pubmed.ncbi.nlm.nih.gov/help/)

### Page navigation
 - [FAQs](#faqs)
- [Search PubMed](#search-pubmed)
- [Display, sort and navigate](#display-sort-navigate)
- [Cite, save, and share](#cite-save-share)
- [Advanced Search](#advanced-search)
- [Other services](#other-services)
- [Appendices](#appendices)
 FAQs Search PubMed Display, sort and navigate Cite, save, and share Advanced Search Other services Appendices
 # PubMed User Guide
 Last update: May 7, 2026
 Visit the [NLM Technical Bulletin](https://techbull.nlm.nih.gov/vivisimo/cgi-bin/query-meta?v:project=technical-bulletin-date&v:sources=technical-bulletin-date&sortby=Date-decending&binning-state=product==PubMed%0a&query=&render.list-show=10) or follow the [PubMed New and Noteworthy RSS feed](https://www.ncbi.nlm.nih.gov/feed/rss.cgi?ChanKey=PubMedNews) for brief announcements highlighting recent enhancements and changes to PubMed.
  ## FAQs
 - How can I [get the full text article](#finding-full-text)? What if the [link to the full text is not working](#broken-linkout-links)?
- How do I [search by author](#author-search)?
- How do I [search by journal name](#journal-search)?
- How do I [find a specific citation](#find-citation)? I have some information such as the author, journal name, and publication year.
- I retrieved too many citations. [How can I focus my search](#narrow-search)?
- I retrieved too few citations. [How can I expand my search](#expand-search)?
- How do I find [consumer health information about a disease or condition?](#consumer-health)
- How do I find [systematic reviews](#systematic-reviews)?
- Are there tools to help with [clinical searches](#clinical-queries) or finding [medical genetics](#medical-genetics) information?
- I’m not finding what I need. [How does a PubMed search work](#automatic-term-mapping)?
- Can you [explain what is shown on the search results](#understanding-docsum)?
- How do I [display an abstract](#display-an-abstract)?
- How can I [save my results](#cite-save-share)?
- Can I [receive email updates when new results are available for my search](#create-an-email-alert)?
- How do I [report an error or duplicate citation in PubMed?](#error-messages)
- How can I [cite an article](#cite) or [export citations to my citation management software program](#citation-management-software)?
- How do I [get a link to bookmark or share my PubMed search](#get-permalink)?
- How can I [download PubMed?](#download-pubmed-data)
- Is there a guide to [NLM resources for MEDLINE/PubMed](https://www.nlm.nih.gov/bsd/pmresources.html)?
- Where can I find [further assistance](#customer-support) and [training](#training)?
  ## Search PubMed
 - [How do I search PubMed?](#how-do-i-search-pubmed)
- [I retrieved too many citations. How can I focus my search?](#narrow-search)
- [I retrieved too few citations. How can I expand my search?](#expand-search)
- [Find a specific citation](#find-citation)
- [Searching by author](#author-search)
- [Searching by journal](#journal-search)
- [Searching by date](#searching-by-date)
- [Filters](#help-filters)
- [Searching for a phrase](#searching-for-a-phrase)
- [Wildcards and truncation](#wildcards)
- [Combining search terms with Boolean operators (AND, OR, NOT)](#combining-with-boolean-operators)
- [Using search field tags](#using-search-field-tags)
- [Proximity searching](#proximity-searching)
 ### How do I search PubMed?
 1. Identify the key concepts for your search.
2. Enter the terms (or key concepts) in the search box.
3. Press the Enter key or click Search.
 For many searches, it is not necessary to use special tags or syntax. PubMed uses multiple tools to help you find relevant results:
 - Best Match sort order uses a state-of-the-art machine learning algorithm to place the most relevant citations at the top of your results.
- An autocomplete feature displays suggestions as you type your search terms. This feature is based on PubMed query log analysis described in " [Finding Query Suggestions for PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815412/)."
- A spell checking feature suggests alternative spellings for search terms that may include misspellings.
- A citation sensor displays suggested results for searches that include terms characteristic of citation searching, e.g., author names, journal titles, publication dates, and article titles.
 ### I retrieved too many citations. How can I focus my search?
 To limit the number of search results:
 - Replace general search terms with more specific ones (e.g., low back pain instead of back pain).
- Include additional terms in your query.
- Use the sidebar filters to restrict results by publication date, full text availability, article type, and more.
 ### I retrieved too few citations. How can I expand my search?
 - On the abstract page for a citation, see the [Similar Articles](#similar-articles) section for a pre-calculated set of additional PubMed citations closely related to that article.
- Remove extraneous or specific terms from the search box.
- Try using alternative terms to describe the concepts you are searching.
 ### Find a specific citation
 Paste the article title into the search box, or enter citation details such as the author, journal name and the year the article was published in the search box and the PubMed citation sensor will automatically analyze your query for citation information to return the correct citation. The citation sensor incorporates a fuzzy matching algorithm and will retrieve the best match even if a search includes an incorrect term. You do not need to use field tags or Boolean operators.
 ### Searching by author
 Enter the author’s last name and initials without punctuation in the search box, and click Search.
 If you only know the author’s last name, use the author search field tag [au], e.g., brody[au].
 Names entered using either the lastname+initials format (e.g., smith ja) or the full name format (john a smith) and no search tag are searched as authors as well as collaborators, if they exist in PubMed.
 Enter a full author name in natural or inverted order, e.g., julia s wong or wong julia s.
 - Prior to 2002, full author names were not included on PubMed citations, so full author name searches will only retrieve citations from 2002 forward, when the full author name was published in the article.
- A comma following the last name for searching is optional. For some names, however, it is necessary to distinguish which name is the last name by using the comma following the last name, e.g., james, ryan.
 Omit periods after initials and put all suffixes at the end, e.g., vollmer charles jr
 Initials and suffixes are not required. If you include a middle initial or suffix, you will only retrieve citations for articles that were published using the middle initial or suffix.
 More information about author searching:
 - To search by author using the search builder, click [Advanced search](https://pubmed.ncbi.nlm.nih.gov/advanced/) and then select Author from the All Fields menu. The author search box includes an autocomplete feature.
- You may click an author link on the abstract display to execute a search for the author in PubMed. Results will display using a ranking algorithm if the author name is computationally similar for additional PubMed citations.
- If an author name includes only [stopwords](#help-stopwords), use the author search field tag [au] to search in combination with other terms, e.g., just by[au] seizure.
- Author names are automatically [truncated](#wildcards) to account for varying initials and designations such as Jr. To turn off the truncation, use double quotes around the author's name with the author search field tag [au], e.g., "smith j"[au].
- Use the search field tag [1au] to search for the first personal author or [lastau] to search for the last personal author name in a citation.
 For additional information on author names in PubMed, please see the journal article, " [Author Name Disambiguation for PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5530597/)."
 ### Searching by journal
 Enter one of the following in the search box:
 - full journal title (e.g., molecular biology of the cell)
- title abbreviation (e.g., mol biol cell)
- ISSN number, a standardized international code (e.g., 1059-1524)
 More information about journal searching:
 - To search by journal using the search builder, click [Advanced search](https://pubmed.ncbi.nlm.nih.gov/advanced/) and then select Journal from the All Fields menu. The journal search box includes an autocomplete feature.
- To find full journal names, use the [NLM Catalog](https://www.ncbi.nlm.nih.gov/nlmcatalog), or mouseover the journal title abbreviation on the citation (available in abstract view). 1. Click Journals in NCBI Databases on the PubMed homepage.
2. Enter the journal name and click Search.
- Use the journal search field tag [ta] to limit your search to the journal only, e.g., gene therapy[ta], scanning[ta]
- Searching with the full journal title or abbreviation is recommended for complete retrieval of indexed items; older citations may not have an ISSN.
- If a journal title or abbreviation includes a special character (e.g., parentheses, brackets, &), enter the title or abbreviation without the special characters. For example, to search by the journal abbreviation j hand surg [am], enter j hand surg am.
- Searching for a journal will automatically map to the official journal title and the title associated with an alternative title, if one exists. To turn off this automatic mapping enter the journal in double quotes and tag with [ta], e.g., "science"[ta].
 A [list of journals](#journal-lists) included in PubMed is available by FTP.
 ### Searching by date
 - [Using the results timeline](#date-search-timeline)
- [Using the search builder](#date-search-builder)
- [Searching by a single date in the search box](#date-search-single)
- [Searching for a date range in the search box](#date-search-range)
- [Searching for a relative date range](#date-search-relative)
 #### Using the results timeline
 Click and drag the sliders on the Publication Date timeline to change the date range for your search.
 Note: The Publication Date timeline counts all publication dates for a citation as supplied by the publisher, e.g., print and electronic publication dates. These dates may span more than one year; for example, an article that was published online in November 2018 and published in a print issue in January 2019. This means the sum of results represented in the timeline may differ from the search results count.
 #### Using the search builder
 1. Click [Advanced search](https://pubmed.ncbi.nlm.nih.gov/advanced/) and use the search builder.
2. Select a date field from the All Fields menu, e.g., Date – Publication, and enter a single date or a date range in the fill-in-the-blank boxes. Month and day are optional. If you want to search for a date range up to the current date, do not edit the ‘Present’ date box.
3. Add the date from the builder to the query box.
4. Once you have finished adding terms to the query box, click Search (or Add to History) to run the search.
 #### Searching by a single date in the search box
 Enter dates using the format yyyy/mm/dd [date field]. The month and day are optional.
 Use a [Boolean operator](#combining-with-boolean-operators) when combining a date with other search terms.
 Example
 Use the Boolean operator AND to limit your search to a specific publication date.
 [cancer AND 2020/06/01[dp]](https://pubmed.ncbi.nlm.nih.gov/?term=cancer+AND+2020%2F06%2F01%5Bdp%5D)
 The available date fields are:
 - [Date of Publication [dp]](#dp) - Includes both print and electronic dates of publication except for electronic publication dates occurring *after* the print publication date.
- Electronic Date of Publication (if applicable) [epdat]
- Print Date of Publication (if applicable) [ppdat]
- [Entry Date [edat]](#edat) - Date used for PubMed processing, such as “Most Recent” sort order.
- [MeSH Date [mhda]](#mhda) - The date the citation was indexed with MeSH terms.
- [Create Date [crdt]](#crdt) - The date the PubMed record was first created.
 #### Searching for a date range in the search box
 Enter date ranges using a colon (:) between each date followed by a [date field].
 Use a [Boolean operator](#combining-with-boolean-operators) when combining a date range with other search terms.
 Example
 Use the Boolean operator AND to limit your search to a date range.
 [heart disease AND 2019/01/01:2019/06/30[dp]](https://pubmed.ncbi.nlm.nih.gov/?term=heart+disease+AND+2019%2F01%2F01%3A2019%2F06%2F30%5Bdp%5D) [influenza AND 2000:2010[dp]](https://pubmed.ncbi.nlm.nih.gov/?term=influenza+AND+2000%3A2010%5Bdp%5D)
 Date range searching includes both print and electronic dates of publication, except for electronic publication dates occurring *after* the print publication date.
 #### Searching for a relative date range
 Use the following format to search for a relative date range:
 - term="last **X** days" **[date field]**
- term="last **X** months" **[date field]**
- term="last **X** years" **[date field]**
 where **X** is the number of days, months or years immediately preceding today’s date and **[date field]** is the date field tag: [dp], [edat] or [crdt].
 The relative date range search for publication dates will also include citations with publication dates after today's date; therefore, citations with publication dates in the future will be included in the results.
 ### Filters
 You can use filters to narrow your search results by [publication date](#filters-publication-date), [text availability](#filters-text-availability), [article attribute](#filters-article-attribute), [article type](#filters-article-type), [article language](#filters-language), [species](#filters-species), [sex](#filters-sex), [age](#filters-age), and [other](#filters-other).
 To apply a filter:
 1. Run a search in PubMed.
2. Click the filter you would like to activate from the sidebar. A check mark will appear next to the activated filter(s).
3. Subsequent searches will be filtered until the selected filters are turned off, or until your browser data is cleared.
4. A green "Filters applied" message will display at the top of the results page when filters are applied to your search.
 The most frequently used filters appear on the sidebar by default. More filters are available for some categories:
 1. Click the "See all … filters" link for a filter category.
2. A pop-up menu will appear showing all available filters for that category.
3. Select the filter(s) you would like to apply to your search. - To cancel your selections, click Cancel or click on the X in the upper right corner to close the pop-up and return to your search results.
4. Click Apply. This will close the pop-up menu, apply the selected filter(s) to your search, and display your selection(s) on the sidebar. - The selected filter(s) will appear on the sidebar in addition to the default filters until your browser cookies are cleared, or you click the "Reset filters display" button at the bottom of the sidebar.
 Clearing applied filters:
 - To turn off all currently applied filters, click the "Clear filters" link in the green "Filters applied" message at the top of the page, or click the "Clear all filters" button at the bottom of the sidebar.
- To turn off individual filters, click the filter on the sidebar.
- To reset the filters shown on the sidebar to the default filter list, click the "Reset filters display" button at the bottom of the sidebar. This will also clear any currently applied filters.
 More information about filters:
 - Many filters under the Additional Filters section limit your search results to citations that have been assigned specific MeSH terms during MEDLINE indexing. This means these filters may exclude some otherwise relevant results that do not include the associated MeSH terms.
- See [Filter search strategies](#filter-strategies) for the equivalent PubMed query for each filter.
- You can activate additional filters with [My Custom filters](https://www.ncbi.nlm.nih.gov/books/NBK53591/). Once activated, these filters will appear at the top of your sidebar under "My Custom Filters".
 #### Publication date filters
 To filter your results by Publication Date, click 1 year, 5 years, 10 years, or enter a custom range. The 1 year, 5 years, and 10 years filters [search for a relative date range](#date-search-relative) based on today's date.
 Publication Date filters include both electronic and print publication dates.
 #### Text availability filters
 To filter your results to only citations that include a link to full text, a link to free full text, or an abstract, click the appropriate selections.
 Alternatively, you may search for citations with links to full text, free full text or include an abstract using the values: full text[sb], free full text[sb], or 'hasabstract'. No search field tag is required for hasabstract. You may also search for all citations with a structured abstract with ‘hasstructuredabstract’.
 Note: Most citations in PubMed to articles published before 1975 do not include abstracts.
 #### Article attribute filter
 Use the associated data filter to find citations that include secondary source databanks and accession numbers. These data links connect to records in other NLM databases (for example: GenBank, ClinicalTrials.gov) or external data repositories (for example: figshare, Dryad). Data links are provided by the publisher and/or by data repositories through participation in NCBI’s LinkOut service.
 #### Article type filters
 Use article type filters to narrow your results based on the content or style of the article, such as: Clinical Trial or Review. Article types are also called publication types.
 You can add more article type filters to the sidebar using the "See all article type filters" link. Currently indexed publication types are available from the filter sidebar; see [Publication Types](#publication-types) for a complete list of all publication types found in PubMed including discontinued publication types, which are not available in the filter sidebar.
 Article type filters typically rely on citations' Publication Type [pt] data. One exception is the [Systematic Review filter](#systematic-reviews), which uses a search strategy in addition to the publication type. See [Filter search strategies: Article type](#filter-strategy-article-type) for more information.
 Publication type data may be supplied by the publisher or assigned during the MEDLINE indexing process. Therefore, these filters may exclude some citations that have not yet completed the MEDLINE indexing process.
 ##### Systematic reviews
 To search for systematic reviews in PubMed, use the Systematic Review article type filter on the sidebar, or enter your search terms followed by AND systematic[sb] in the search box. For example, lyme disease AND systematic[sb].
 The Systematic Review filter uses a [search strategy](https://www.nlm.nih.gov/bsd/pubmed_subsets/sysreviews_strategy.html) in addition to the Systematic Review publication type [pt] to find systematic reviews in PubMed. To limit your search to only those citations with the Systematic Review publication type, use the [publication type](#publication-types) search tag[pt], i.e., systematic review[pt]; however, this may exclude some relevant citations that have not yet completed the MEDLINE indexing process.
 #### Article language filters
 Language filters restrict your search to articles published in the selected language(s). You can add more language filters to the sidebar using the "See all article language filters" link.
 By default, PubMed displays English language titles and abstracts when provided by the publisher. Check the Abstract display for links to view the abstract in other languages (when available).
 Article language filters are located under the Additional Filters section.
 #### Species filters
 Species filters restrict your results to human or animal studies:
 - Humans
- Other Animals
 Species filters limit your search results to citations that were assigned [the associated MeSH terms for species](#filter-strategy-species) during MEDLINE indexing. These filters will typically deliver precise results but can also exclude some otherwise relevant results that do not include the associated MeSH terms on their PubMed record. For example, some citations are not indexed for MEDLINE, such as preprints, Online Ahead of Print citations, and U.S. government-funded research published in journals that are not indexed for MEDLINE.
 Species filters are located under the Additional Filters section.
 #### Sex filters
 Sex filters restrict your search results to a specific sex for an animal or human study:
 - Female
- Male
 Sex filters limit your search results to citations that were assigned [the associated MeSH terms for sex](#filter-strategy-sex) during MEDLINE indexing. These filters will typically deliver precise results but can also exclude some otherwise relevant results that do not include the associated MeSH terms on their PubMed record. For example, some citations are not indexed for MEDLINE, such as preprints, Online Ahead of Print citations, and U.S. government-funded research published in journals that are not indexed for MEDLINE.
 Sex filters are located under the Additional Filters section.
 #### Age filters
 Age filters restrict results to a specific age group for a human study.
 Age filters limit your search results to citations that were assigned [the associated MeSH terms for age](#filter-strategy-age) during MEDLINE indexing. These filters will typically deliver precise results but can also exclude some otherwise relevant results that do not include the associated MeSH terms on their PubMed record. For example, some citations are not indexed for MEDLINE, such as preprints, Online Ahead of Print citations, and U.S. government-funded research published in journals that are not indexed for MEDLINE.
 Age filters are located under the Additional Filters section. You can add more age filters to the sidebar using the "See all age filters" link.
 Age filters include:
 - Child: birth-18 years
- Newborn: birth-1 month
- Infant: birth-23 months
- Infant: 1-23 months
- Preschool Child: 2-5 years
- Child: 6-12 years
- Adolescent: 13-18 years
- Adult: 19+ years
- Young Adult: 19-24 years
- Adult: 19-44 years
- Middle Aged + Aged: 45+ years
- Middle Aged: 45-64 years
- Aged: 65+ years
- 80 and over: 80+ years
 #### Other filters & more subsets
 ##### Exclude preprints
 The Exclude preprints filter will exclude preprint citations from your search results. Alternatively, you can exclude preprints from your search results by including NOT preprint[pt] at the end of your query.
 See [Preprints](#coverage-preprints) for more information about preprint citations in PubMed.
 The Exclude preprints filter is located under the Additional Filters section.
 ##### MEDLINE Subset
 The MEDLINE filter limits results to citations that are indexed for MEDLINE. To use this filter in a query, add medline[sb] to your search.
 The MEDLINE filter is located under the Additional Filters section.
 ##### PubMed Central Subset
 To restrict retrieval to citations that have a free full text article available in [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/) (PMC), search "pubmed pmc"[sb].
 Use the [PMID/PMCID/NIHMSID Converter](https://www.ncbi.nlm.nih.gov/pmc/pmctopmid/) to convert IDs for publications referenced in PubMed and PMC. To retrieve citations that include an NIHMS ID use the query, hasnihmsid.
 ##### Citation Status Subsets
 The citation status indicates the internal processing stage of an article in the PubMed database (see [PubMed Citation Status Subsets](#citation-status-subsets)).
 To search for a particular citation status, enter one of the search terms below followed by the [sb] search tag:
 - publisher
- inprocess
- medline
- pubmednotmedline
 To search for the total number of PubMed citations, enter all[sb] in the search box.
 ##### Ahead of Print Citations
 Publishers may submit citations for articles that appear on the web prior to their publication in final or print format. To search for these ahead-of-print citations, enter pubstatusaheadofprint.
 ### Searching for a phrase
 Many phrases are recognized by the subject translation table used in [PubMed's Automatic Term Mapping (ATM)](#automatic-term-mapping). For example, if you enter fever of unknown origin, PubMed recognizes this phrase as a MeSH Term.
 You can bypass ATM and search for a specific phrase using the following formats:
 - Enclose the phrase in double quotes: "kidney allograft" - If you use quotes and the phrase is not found in the [phrase index](#phrase-index), the quotes are ignored and the terms are processed using automatic term mapping. The message "Quoted phrase not found in phrase index" will display at the top of your search results.
- Use a search tag: kidney allograft[tw] - If you use a search tag and the phrase is not found in the [phrase index](#phrase-index), the phrase will be broken into separate terms. For example, "psittacine flight" is not in the phrase index, so a search for psittacine flight[tw] is broken up and translated as: ((("psittaciformes"[MeSH Terms] OR "psittaciformes"[All Fields]) OR "psittacine"[All Fields]) OR "psittacines"[All Fields]) AND "flight"[Text Word]
- Terms tagged with [All] are processed using Automatic Term Mapping (ATM). This means phrases tagged with [All] matching a MeSH term or entry term will also include those terms in the search translation. If you want to prevent such mapping, put the phrase inside double quotes: "heart attack"[All]
- Use a hyphen: kidney-allograft - If you use a hyphen and the phrase is not found in the [phrase index](#phrase-index), the search will not return any results for that phrase.
- Hyphenated phrases matching a MeSH term or entry term will include those terms in the search translation. If you want to prevent such mapping, put the hyphenated phrase inside double quotes: "heart-attack"
 When you enter search terms as a phrase, PubMed will not perform automatic term mapping that includes the MeSH term and any specific terms indented under that term in the MeSH hierarchy. For example, "health planning" will include citations that are indexed to the MeSH term, Health Planning, but will not include the more specific terms, e.g., Health Care Rationing, Health Care Reform, Health Plan Implementation, that are included in the automatic MeSH mapping.
 #### Phrase index
 PubMed uses a phrase index to provide phrase searching. To browse the phrase index, use the [Show Index feature](#browsing-show-index) included in the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) builder: select a search field, enter the beginning of a phrase, and then click Show Index.
 ##### Quoted phrase not found
 Phrases may appear in a PubMed record but not be in the phrase index. To search for a phrase that is not found in the phrase index, use a [proximity search](#proximity-searching) with a distance of 0 (e.g., ["cognitive impairment in multiple sclerosis"[tiab:~0]](https://pubmed.ncbi.nlm.nih.gov/?term=%22cognitive+impairment+in+multiple+sclerosis%22%5Btiab%3A%7E0%5D)); this will search for the quoted terms appearing next to each other, in any order.
 Automated processes regularly add new phrases to the index based on standard criteria such as phrase frequency and length. If you would like to request a phrase be added to the phrase index, please [write to the NLM Help Desk](https://support.nlm.nih.gov/support/create-case/).
 ### Wildcards and truncation
 To create a wildcard search in PubMed, use an asterisk (*) to substitute for 0 or more characters in a term or phrase. For example, ["vaccin* schedul*"](https://pubmed.ncbi.nlm.nih.gov/?term=%22vaccin*+schedul*%22) will find results that include "vaccine schedule", "vaccine scheduling", "vaccination schedules", and more variations of this phrase.
 - A wildcard (*) substitutes for 0 or more characters
- Terms must begin with at least 4 characters before the first wildcard: colo*
- Multiple wildcards can be used in the same term: organi*ation*
 To use wildcards in a phrase search, use any of the following formats:
 - Enclose the phrase in double quotes: "breast feed*"
- Use a search tag: breast feed*[tiab]
- Use a hyphen: breast-feed*
 Phrases can include more than one wildcard: ["colo* cancer*"](https://pubmed.ncbi.nlm.nih.gov/?term=%22colo*+cancer*%22)
 Wildcards turn off [Automatic Term Mapping (ATM)](#automatic-term-mapping) and the process that includes the MeSH term and any specific terms indented under that term in the MeSH hierarchy. For example, "heart attack*" will not map to the MeSH term Myocardial Infarction or include any of the more specific terms, e.g., Myocardial Stunning; Shock, Cardiogenic.
 Wildcards can be particularly useful to capture variations of a phrase, since phrase searches already turn off automatic term mapping and do not otherwise include alternate spellings or singular/plural forms. For example, searching for "tumo*r associated macrophage*" retrieves variations of the phrase such as: "tumor associated macrophage", "tumour associated macrophage", "tumor associated macrophages", and "tumour associated macrophages".
 ### Combining search terms with Boolean operators (AND, OR, NOT)
 PubMed applies an AND operator between concepts, e.g., "vitamin c common cold" is translated as vitamin c AND common cold. Enter Boolean operators in uppercase characters to combine or exclude search terms:
 - AND retrieves results that include all the search terms.
- OR retrieves results that include at least one of the search terms.
- NOT excludes the retrieval of terms from your search.
 PubMed processes searches in a left-to-right sequence. Use parentheses to "nest" concepts that should be processed as a unit and then incorporated into the overall search.
 - PubMed uses automatic term mapping to identify concepts. For example, for the search air bladder fistula, PubMed will search "air bladder" as a phrase. If you do not want this automatic phrase parsing, enter each term separated by the Boolean operator AND, e.g., air AND bladder AND fistula.
- [Search Details](#viewing-search-details) show how a search was translated.
 ### Using search field tags
 You can search for a term in a specific field by including a [search field tag](#search-tags) after the term; for example, UCLA[ad] will search for the term “UCLA” in the affiliation field only.
 More information about using search field tags:
 - The search field tag must be enclosed in square brackets.
- Case and spacing do not matter: crabs [mh] = Crabs[mh].
- Search field tags turn off [Automatic Term Mapping (ATM)](#automatic-term-mapping), limiting your search to the specified term only.
- Using a search field tag after multiple terms will attempt to [search those terms as a phrase](#searching-for-a-phrase): kidney allograft[tiab].
- To search multiple terms in the same field, each term must be tagged individually: covid-19[ti] vaccine[ti] children[ti].
- The [Advanced Search builder](#advanced-search) can help you search for terms in specific fields and build large, complex search strings.
 #### Search field tags
 - [Affiliation [ad]](#ad)
- [All Fields [all]](#all)
- [Article Identifier [aid]](#aid)
- [Author [au]](#au)
- [Author Identifier [auid]](#auid)
- [Book [book]](#book)
- [Comment Correction Type](#comment-correction)
- [Completion Date [dcom]](#dcom)
- [Conflict of Interest Statement [cois]](#cois)
- [Corporate Author [cn]](#cn)
- [Create Date [crdt]](#crdt)
- [EC/RN Number [rn]](#rn)
- [Editor [ed]](#ed)
- [Entry Date [edat]](#edat)
- [Filter [filter] [sb]](#filter-sb)
- [First Author Name [1au]](#1au)
- [Full Author Name [fau]](#fau)
- [Full Investigator Name [fir]](#fir)
- [Grants and Funding [gr]](#gr)
- [Investigator [ir]](#ir)
- [ISBN [isbn]](#isbn)
- [Issue [ip]](#ip)
- [Journal [ta]](#ta)
- [Language [la]](#la)
- [Last Author Name [lastau]](#lastau)
- [Location ID [lid]](#lid)
- [MeSH Date [mhda]](#mhda)
- [MeSH Major Topic [majr]](#majr)
- [MeSH Subheadings [sh]](#sh)
- [MeSH Terms [mh]](#mh)
- [Modification Date [lr]](#lr)
- [NLM Unique ID [jid]](#jid)
- [Other Term [ot]](#ot)
- [Owner](#owner)
- [Pagination [pg]](#pg)
- [Personal Name as Subject [ps]](#ps)
- [Pharmacological Action [pa]](#pa)
- [Place of Publication [pl]](#pl)
- [PMCID and MID](#pmcid-mid)
- [PMID [pmid]](#pmid)
- [Publication Date [dp]](#dp)
- [Publication Type [pt]](#pt)
- [Publisher [pubn]](#pubn)
- [Secondary Source ID [si]](#si)
- [Subset [sb]](#sb)
- [Supplementary Concept [nm]](#nm)
- [Text Words [tw]](#tw)
- [Title [ti]](#ti)
- [Title/Abstract [tiab]](#tiab)
- [Transliterated Title [tt]](#tt)
- [Volume [vi]](#vi)
 ### Proximity searching
 You can use proximity searching to search for multiple terms appearing in any order within a specified distance of one another in the [Title], [Title/Abstract], or [Affiliation] fields.
 To create a proximity search in PubMed, enter your terms using the following format:
 **"search terms"[field:~N]**
 - **Search terms** = Two or more words enclosed in double quotes. - There is no limit to the number of words you can search together in proximity; however, the more terms you enter, the more restrictive your search becomes. Using the [Boolean operator](#combining-with-boolean-operators) AND to combine terms may be more appropriate than combining many terms into one proximity search.
- **Field** = The search field tag for [Title], [Title/Abstract], or [Affiliation] fields.
- Proximity searching is only available in the [Title](#ti), [Title/Abstract](#tiab), and [Affiliation](#ad) search fields.
- You can use the full search field tags [Title], [Title/Abstract], and [Affiliation], or the abbreviated versions [ti], [tiab], and [ad].
- **N** = The maximum number of words appearing between your search terms. - What N value to use will depend on your search. Try changing the N value and comparing the results to find what works best for your search.
- A higher N creates a broader, more comprehensive search; this will typically retrieve more results overall, but some of these results may be less relevant. Using the [Boolean operator](#combining-with-boolean-operators) AND to combine terms may be more appropriate than proximity searching with a large N value.
- A lower N creates a narrower, more precise search; this will typically retrieve fewer results that are highly relevant, but may exclude other relevant results.
- If N=0, the quoted terms will appear next to each other--with no other words in between.
- For the affiliation field only, an N value of 1,000 or less will search for the double quoted terms together within the same affiliation, rather than spread across all affiliations on the record. See [Affiliation [ad] for an example proximity search in the affiliation field](#ad) and more information about searching for affiliations.
 More information about proximity searching:
 - Results will include your quoted terms in any order. If you would like to search for an exact phrase with terms appearing in a specific order, use a phrase search instead.
- Automatic Term Mapping is not applied to the quoted terms.
- Proximity searching is not compatible with [wildcards](#wildcards) (*). If the double quoted terms in a proximity search include a wildcard (*), the proximity operator will be ignored.
- You can combine proximity searches with other terms using Boolean operators; for example, "hip pain"[Title:~4] AND stretching
- Booleans and [stopwords](#help-stopwords) included in quoted terms for proximity search are searched like regular keywords.
 Example
 Search PubMed for citations with the terms "rationing" and "healthcare" appearing within 2 words of each other--in any order--in the Title field:
 ["rationing healthcare"[Title:~2]](https://pubmed.ncbi.nlm.nih.gov/?term=%22rationing+healthcare%22%5BTitle%3A%7E2%5D) Search results may include: rationing healthcare, healthcare rationing, rationing of healthcare, rationing in healthcare, rationing universal healthcare, rationing strategies in healthcare, rationing limited healthcare… and more.
 Example
 Search PubMed for citations with the terms "patient," "physician," and "relationship" appearing next to each other—in any order—in the Title/Abstract fields:
 ["patient physician relationship"[Title/Abstract:~0]](https://pubmed.ncbi.nlm.nih.gov/?term=%22patient+physician+relationship%22%5BTitle%2FAbstract%3A%7E0%5D) Since N=0, the quoted terms must appear next to each other with no other words in between them, although they can still appear in any order.
 Example
 You can build queries that combine proximity searches with other terms using Boolean operators (AND, OR, NOT):
 ["accidental overdose"[Title/Abstract:~3] AND acetaminophen](https://pubmed.ncbi.nlm.nih.gov/?term=%22accidental+overdose%22%5BTitle%2FAbstract%3A%7E3%5D+AND+acetaminophen) ["patient physician relationship"[tiab:~0] OR “patient doctor relationship”[tiab:~0]](https://pubmed.ncbi.nlm.nih.gov/?term=%22patient+physician+relationship%22%5Btiab%3A%7E0%5D+OR+%22patient+doctor+relationship%22%5Btiab%3A%7E0%5D)
  ## Display, Sort, and Navigate
 - [Understanding your search results](#understanding-docsum)
- [Display an abstract](#display-an-abstract)
- [Changing the display format of search results](#changing-display-format)
- [Showing more results](#show-more-results)
- [Sorting your results](#sort-results)
- [Finding the full text article](#finding-full-text)
- [Figures](#figures)
- [Similar articles](#similar-articles)
- [Cited by](#cited-by)
- [References](#help-references)
- [Grants and funding](#grants-and-funding)
- [Navigating searches with more than 10,000 results](#10k-results)
- [Discovering related data in NCBI databases](#related-data)
- [Find related resources using LinkOut](#using-linkout)
- [Reporting broken or problem links](#broken-linkout-links)
 ### Understanding your search results
 Citations are initially displayed 10 items per page and sorted by Best Match.
 By default, PubMed search results are displayed in a summary format and include snippets from the citation abstract. Snippets and highlighted terms are selected based on relatedness to your query.
 To see the abstract for an individual citation, click the title of the citation to go to its abstract page.
 Journal names are shown using the journal title abbreviation. When viewing citations in Abstract format, you can mouseover a journal’s title abbreviation to display the full journal name.
 ### Display an abstract
 Click the title of the citation to go to its abstract page, or change the search results display to Abstract format using the Display options button in the upper right corner of the search results page.
 PubMed may include non-English abstracts if supplied by the publisher. The abstract text defaults to English when a citation has an accompanying non-English abstract. Links to display the additional language(s) are available on the Abstract display. To retrieve citations with non-English abstracts, use the query hasnonenglishabstract.
 ### Changing the display format of search results
 Results are displayed in the summary format by default, except a single citation result will go directly to the abstract page. You can change the results format using the Display options button:
 1. Click the Display options button in the upper right corner of the search results page
2. Select the display format you would like to use
3. Results will be displayed in the new format
 Selecting one or more items and changing the display format will display only the selected result(s) in the new format.
 By default, the summary format includes snippets from the citation abstract. You can turn off snippets under Display options by deselecting Abstract snippets.
 Additional display preferences are available when logged into My NCBI, such as choosing a highlighting color or automatically displaying search results in Abstract format. Please see our online help for more information:
 - [My NCBI Help: Using Preferences](https://www.ncbi.nlm.nih.gov/books/NBK53593/)
- [My NCBI Help: What is My NCBI?](https://www.ncbi.nlm.nih.gov/books/NBK3842/)
 ### Showing more results
 The results page indicates the total number of items retrieved.
 Ten items are displayed per page by default. You can change the number of items displayed per page using the Display options button:
 1. Click the Display options button in the upper right corner of the search results page
2. Select the number of items to display per page: 10, 20, 50, 100, or 200
3. Your selection will be active for subsequent searches until your browser cookies are cleared.
 Click "Show more" to display the next page of results, or click "Jump to page" to navigate directly to a specific page of results.
 ### Sorting your results
 The default sort order in PubMed is Best Match. You can use the "Sort by" drop-down menu at the top of the search results page to change the sort order.
 If you change the sort order, your new selection will be active for subsequent searches until your browser cookies are cleared.
 #### Sort orders
 You can sort your search results by:
 - Best Match: The Best Match sort order is based on an algorithm that analyzes each PubMed citation found with your search terms. For each search query, "weight" is calculated for citations depending on how many search terms are found and in which fields they are found. In addition, recently-published articles are given a somewhat higher weight for sorting. The top articles returned by the weighted term frequency algorithm above are then re-ranked for better relevance by a new machine-learning algorithm. Please see the [Algorithm for finding best matching citations in PubMed](#best-match-algorithm) for more information.
- Most Recent: Citations sorted by Most Recent are displayed in reverse date added order: last in, first out. The Most Recent date is the date a record was initially added to PubMed, not the publication date. The secondary sort is PMID.
- Publication Date: Citations sorted by Publication Date are displayed in reverse chronological order: newest to oldest. Citations with more than one publication date, such as electronic and print, are sorted by their earliest publication date. Publication dates without a month are set to January, multiple months (e.g., Oct-Dec) are set to the first month, and dates without a day are set to the first day of the month. Dates with a season are set as: winter = January, spring = April, summer = July and fall = October.
- First Author: Citations are sorted alphabetically by first author name. The secondary sort order within a group with the same first author is PMID.
- Journal: Citations are sorted alphabetically by journal name. The secondary sort order within a group with the same journal name is PMID.
 #### Reverse sort order
 - When sorting by Most Recent, Publication Date, First Author, or Journal, you can reverse the sort order by clicking the up/down arrow next to the selected sort option to toggle between ascending or descending order.
- The reverse sort option will not display when Best Match sort order is selected.
 #### Computed author sort
 Clicking an author name link on the abstract display runs a search for the author in PubMed. If an author name is computationally similar with an author name for additional PubMed citations, the results will display those citations first, in ranked order, followed by the non-similar citations. Author name disambiguation details are available in [Liu W and Wilbur WJ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5530597/).
 ### Finding the full text article
 PubMed records contain citation information (e.g., title, authors, journal, publication date) and abstracts of published articles and books. PubMed search results do not include the full text of the journal article, but the abstract view in PubMed includes links to the full text from other sources when available, such as the publisher’s website or the [PubMed Central (PMC)](https://www.ncbi.nlm.nih.gov/pmc/) database. The full text journal site may require a fee or subscription, however online journals sometimes provide free access. Access may also be available through your organization, or local medical library.
 You may be able to obtain free copies of full text articles in these ways:
 #### Free full text filter
 On the filter sidebar, click "Free full text" to narrow results to resources that are available for free on the web, including PubMed Central, Bookshelf, and publishers' websites. Alternately, include free full text[Filter] in your query.
 #### PubMed Central
 When full text is available in [PubMed Central (PMC)](https://www.ncbi.nlm.nih.gov/pmc/), the "Free in PMC" icon will appear on the citation's abstract display under Full Text Links. Click the icon to view the article in PMC.
 PubMed Central (PMC) is the U.S. National Institutes of Health (NIH) free digital archive of biomedical and life sciences journal literature.
 #### From the publisher
 Journal publishers or related organizations may provide access to articles for free, for free after registering as an individual or guest, or for a fee. When provided by the publisher or other organization, icons linking to these sources can be found on the citation's abstract display under the "Full Text Links" and/or "LinkOut" sections. Icons will often indicate free full text when the article is available for free.
 Note: When you click a full text icon or link in PubMed, you leave PubMed and are directed to the full text at an external provider's site. NCBI does not hold the copyright to this material, and cannot give permission for its use. Users should review all copyright restrictions set forth by the full text provider before reproducing, redistributing, or making commercial use of material accessed through LinkOut.
 Please see the [Copyright and Disclaimers](https://www.ncbi.nlm.nih.gov/About/disclaimer.html) page for additional information.
 #### If you are affiliated with a hospital, university, or other institution
 Your local medical library is your best option. If you see icons for your library on the abstract view this indicates that your library provides a link to the article, has the journal in its collection, or may otherwise obtain the article for you through interlibrary loan. If your library does not have access to the article you need, ask a librarian about ordering the article from another institution.
 #### Local library
 Some local libraries have copies of medical journals or can get a copy of an article for you. Ask your local librarian about inter-library loan options and fees.
 ### Figures
 PubMed abstracts include figures when the full text article is available in [PubMed Central (PMC)](#coverage-pmc). Click the thumbnail to view a larger version of the image, caption, and link to the figure and copyright information in PMC. Figures correspond to the article version archived in PMC, which in some cases may be the [Author Manuscript](https://pmc.ncbi.nlm.nih.gov/about/authorms/) version.
 ### Similar articles
 The abstract page for a citation includes links to PubMed citations for similar articles. The "See all similar articles" link will retrieve a pre-calculated set of PubMed citations that are closely related to the selected article:
 - Similar articles are displayed in ranked order from most to least relevant, with the "linked from" citation displayed first.
- Similar articles are generated by comparing words from the title, abstract, and MeSH terms using a word-weighted algorithm.
- Filters are not activated for similar articles.
- You can refine the list of similar articles using your search [History](#combining-history), where the similar articles retrieval is represented as a list of PMIDs. Use this search number in a search. Refining the list removes the ranked order and may remove citations that are most relevant.
 See [Computation of similar articles](#computation-of-similar-articles) for more information.
 ### Cited by
 PubMed abstracts include links to other publications citing the current item. For items with more than 5 citing publications linked in PubMed, a "See all 'cited by' articles" link will retrieve the list of citing publications with the cited article displayed as the first result.
 "Cited by" is generated using data submitted by publishers and from NCBI resources, when available. "Cited by" may not be a complete list of works citing a particular item.
 ### References
 PubMed abstracts include references when available. Reference lists are available for citations to full text articles included in the open access subset of PMC and for citations where the publisher supplied references in the citation data sent to PubMed.
 ### Grants and funding
 PubMed displays grant numbers, contract numbers, and intramural research identifiers that have been associated with a publication by:
 1. Publishers when depositing data in PubMed and PubMed Central;
2. An author, principal investigator, or project director when: - depositing a manuscript through the [NIH Manuscript Submission (NIHMS)](https://www.nihms.nih.gov/) or [Europe PMC Plus](https://plus.europepmc.org/) system; or
- when adding a publication to [My Bibliography](https://www.ncbi.nlm.nih.gov/books/NBK53595/#mybibliography.Associating_Funding_to_yo); and/or
3. NLM text mining and indexing processes.
 A grant award or contract may be acknowledged in an article and, therefore, displayed in PubMed, for various reasons, including support for activities that contributed directly to the publication as well as support for the generation of an underlying dataset or another shared resource. Additionally, some articles may not explicitly acknowledge intramural research support, yet the authors may be affiliated with a funding agency and may have associated their intramural support with a PubMed record at the time of manuscript deposit to PMC.
 Funding information in PubMed is collected in or converted to a [standardized format](#gr) when possible to enable broad discovery and impact monitoring. For example, if a publication acknowledges support from NIH grant number 1R01 GM987654-01-A1 or GM987654 or ROI GM987654 in a publication, in PubMed the funding information would be normalized to R01 GM987654, consistent with [NIH requirements](https://grants.nih.gov/policy-and-compliance/policy-topics/federal-funding#proper-grant-number-format) for proper grant number format. Funding associations made in a manuscript submission, grant reporting, or indexing system use standardized project identifiers provided to NLM by the organization administering the funding. To learn about searching funding information, see the search field section on [Grants and funding [gr]](#gr).
 The scope of funding information included in PubMed has expanded over time to support the public access policies of [NIH and other funding organizations](https://www.ncbi.nlm.nih.gov/pmc/about/public-access/). Since 1981, NLM has included grant or contract numbers or both that designate financial support by any agency of the United States Public Health Service (PHS), including NIH. Until 2000, only up to three grant numbers were included. Beginning in March 2006, funding information was expanded in PubMed to include grant, contract, and intramural funding assertions made in NIHMS and My Bibliography to support the NIH Public Access Policy. Publishers have been able to supply funding information directly to PubMed since January 2017.
 #### Reporting funding information errors
 Some publications may be inadvertently linked to the wrong funding information. For example, the association of a publication to NIH-funded extramural research requires that the author(s) acknowledge NIH support in the article and that the acknowledgement be in a form that can be readily associated with a specific grant or contract. Variations in the format used to cite NIH funding may lead to either an inability to make an association or erroneous matches of publications to grants and contracts.
 If you identify an error in funding information associated with a PubMed record, please contact the [NLM help desk](https://support.nlm.nih.gov/support/create-case/). NLM will not remove funding associations that reflect the acknowledged funding in the article without a published correction to ensure alignment with the scientific record. If an award association was provided by the author, principal investigator, or project director in My Bibliography or the NIHMS for formal NIH progress and public access compliance reporting, removing the association requires the principal investigator be notified and confirm the lack of direct support.
 ### Navigating searches with more than 10,000 results
 PubMed can display up to 10,000 results. The following options can help you navigate searches with more than 10,000 results:
 - [Reverse the sort order](#reverse-sort) to see the last results first.
- Divide the result set into smaller chunks using the [results timeline](#date-search-timeline) or [custom date range filter](#filters-publication-date).
- [Adjust your search](#narrow-search) to retrieve fewer results.
- For programmatic use and bulk downloads, [PubMed data is available via FTP](#download-pubmed-data).
 ### Discovering related data in NCBI databases
 When available, links to other related NCBI databases are included on a citation's Abstract page under the Related information section. The complete list of database options is provided in [Entrez Link Descriptions](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html#pubmed).
 MEDLINE indexed citations include additional supplemental information on the Abstract page such as MeSH terms, publication types, and substances with links to search for these data in PubMed and the MeSH Database.
 To simultaneously search all NCBI databases, use the [NCBI Search page](https://www.ncbi.nlm.nih.gov/search/).
 ### Find related resources using LinkOut
 Most PubMed records include [LinkOut](https://www.ncbi.nlm.nih.gov/projects/linkout/) resources to a variety of websites including publishers, aggregators, libraries, biological databases, and sequence centers. LinkOut resources link to providers’ sites to obtain the full text of articles or related information, e.g., consumer health. There may be a charge to access the text or information from a provider's site.
 To view LinkOut resources, navigate to the LinkOut section at the end of an individual citation's abstract page.
 To find citations with links to free full text articles, apply the "Free full text" filter to your search results.
 To find citations with links to full text articles, enter search terms followed by AND full text[sb].
 More information about Links:
 - LinkOut resource categories such as "free full text" have been selected by the LinkOut provider.
- The current list of [LinkOut providers](https://www.ncbi.nlm.nih.gov/projects/linkout/journals/htmllists.cgi?type_id=16) is available.
- A publisher's icon link may display on the abstract format if they have electronically provided their citation data to NCBI. Links are only available for publishers that are participating in LinkOut; publishers are responsible for providing working links.
 ### Reporting broken or problem links
 LinkOut links are supplied by the LinkOut providers. Publishers who electronically supply their data to PubMed may include an icon that links to a site providing the full text. Corrections and changes to links are made by the providers and are their responsibility.
 To report problem links or inquire about online journal subscriptions, contact the provider directly. Contact information is typically available at a provider's web site.
  ## Cite, Save, and Share
 - [Save citations temporarily using the Clipboard](#save-citations-to-clipboard)
- [Save citations indefinitely using My NCBI Collections](#save-citations-to-collection)
- [Save citations as a text file](#save-citations-to-file)
- [Cite an article](#cite)
- [Export citations into citation management software](#citation-management-software)
- [Email citations](#email-citations)
- [Create an email alert for a search](#create-an-email-alert)
- [Create an RSS feed for a search](#create-an-rss-feed)
- [Print your search results](#print-search-results)
- [Get a permalink to bookmark or share your search](#get-permalink)
- [Download PubMed data](#download-pubmed-data)
 ### Save citations temporarily using the Clipboard
 The Clipboard provides a place to collect up to 500 items from one or more searches. Items saved to the Clipboard are stored in your browser cookies and will expire after 8 hours of inactivity. If you would like to save items for longer than 8 hours or to view on another device, please use [Send to: Collections](#save-citations-to-collection).
 To add items to the Clipboard:
 1. Use the check boxes to select items from your search results. To save all results (up to a maximum of 500), do not tick any check boxes.
2. Use the Send to button and choose Clipboard.
3. Selected items will be added to the Clipboard. - If no items were selected, a drop-down menu of options will display where you may add selected items, all results on the page, or all results (up to a maximum limit of 500 citations) to the Clipboard.
- An individual item can also be added to the Clipboard from its abstract page.
4. To view your selections, click the Clipboard link under the Search bar. This link will only appear after one or more items have been added to the Clipboard; the link is not present when the Clipboard is empty.
 To delete items from the Clipboard:
 - On the Clipboard page, click "Remove from Clipboard" below each item to delete the item from the Clipboard.
- Select one or more items using the check boxes next to each item, then click "Remove selected items."
- To delete all items from the Clipboard, click "Remove all."
 More information about the Clipboard:
 - Citations added to the Clipboard are marked with the message "Item in Clipboard" in search results.
- The maximum number of items that can be sent to the Clipboard is 500. If you select Clipboard from send to without selecting citations, PubMed will add all (up to 500 citations) of your search results to the Clipboard.
- The Clipboard will not add a citation that is currently in the Clipboard; it will not create duplicate entries.
- Your web browser must accept [cookies](#help-cookies) to use the Clipboard.
- Citations in the Clipboard are represented by the search number #0, which may be used in Boolean search statements. For example, to limit the citations you have collected in the Clipboard to English language articles, use the following search: #0 AND english [la]. This does not affect or replace the Clipboard contents.
 ### Save citations indefinitely using My NCBI Collections
 Search results can be saved in My NCBI using the [Collections](https://www.ncbi.nlm.nih.gov/books/NBK53590/) feature. There is no limit to the number of collections you may store in My NCBI. In addition, collections can be made public to share with others.
 To save results to a new collection:
 1. Sign into My NCBI.
2. Run a search in PubMed.
3. Use the check boxes to select items from your search results or Clipboard. To save all results (up to a maximum of 1,000), do not tick any check boxes.
4. Use the Send to button and choose Collections.
5. Selected items will be added to a Collection. - If no items were selected, a drop-down menu of options will display where you may add selected items, all results on the page, or all results (up to a maximum limit of 1,000 citations) to a Collection.
- An individual item can also be added to a Collection from its abstract page.
6. Choose Create a new collection.
7. Name your collection using a short, meaningful title. The name must be unique and less than 100 characters. Identical names for different Collections are not allowed.
8. Click Add to finish.
 As you continue to build collections, you may want to add new items to an existing collection. To add search results to an existing collection:
 1. Follow steps 1 - 5 above. Add to an existing collection will be the default selection.
2. Use the pull-down menu to choose a collection.
3. Click Add to finish.
 For more information on viewing, sorting, editing, merging, sharing, and deleting collections, see [Collections](https://www.ncbi.nlm.nih.gov/books/NBK53590/) in My NCBI Help.
 ### Save citations as a text file
 Use the Save button to download citations to a text file.
 1. Use the check boxes to select citations from your search results or Clipboard. You may move to other pages to continue your selections. If you do not make any selections, you can choose to save “All results on this page” or “All results” from the Save menu.
2. Click Save and choose a Selection and Format from the menu that appears. - Selection: The citations you would like to save. - Selection: The number of selected items will be shown, for example: Selection (87).
- All results on this page
- All results (up to a maximum of 10,000 citations)
- Format: Summary (text), [PubMed](#pubmed-format), PMID list, Abstract (text), or CSV
3. Click Create file.
4. Your web browser will prompt you to save the file on your computer.
 More information about saving citations to a file:
 - Saving a large set of results may take several minutes.
- To save citations in HTML format, use the "Save" or "Save as" function of your browser and change the file extension to html. When saving as HTML, only those citations displayed on the page will be saved; therefore, consider [showing more results](#show-more-results).
 ### Cite an article
 The Cite button makes it easy to retrieve styled citations that you can copy and paste into a document, or download an .nbib file to use with your reference manager software.
 Using the Cite button for an item will open a pop-up window where you can copy the citation formatted in four popular styles: AMA (American Medical Association), MLA (Modern Language Association), APA (American Psychological Association), or NLM (National Library of Medicine). You can also download the citation as an .nbib file, which most bibliographic reference management software can import.
 Note: In all citation styles, there are certain capitalization rules that machines cannot handle. For example, there is no way to identify proper nouns, acronyms, abbreviations, etc., that is 100% accurate and complies with all rules at all times. Capitalization of article titles and other citation elements should be checked for compliance with a particular reference style when required.
 To export multiple citations: follow the instructions for [saving citations as a text file](#save-citations-to-file) and choose the format Summary (text) to save a list of citations in NLM style, or follow the instructions to [export citations into your citation management software program](#citation-management-software).
 ### Export citations into citation management software
 Use Send to: Citation Manager to export citations as an .nbib file that can be used by many citation management programs:
 1. Use the check boxes to select citations from your search results or Clipboard. You may move to other pages to continue your selections. Alternately, you can choose to save all results on this page or all results from the Send to: Citation Manager menu.
2. Click Send to and choose Citation Manager.
3. Confirm the citations you want to export: selection, all results on this page, or all results (up to a maximum of 10,000).
4. Click Create file.
5. Your web browser will prompt you to save the file on your computer.
6. Import this saved file into your citation management program.
 You can also download an .nbib file for individual citations using the Cite button.
 Questions regarding citation management software should be directed to the respective companies.
 ### Email citations
 [Sign in to My NCBI](https://www.ncbi.nlm.nih.gov/myncbi/) to use this feature. PubMed emails will be sent to the email address associated with the My NCBI account.
 1. Use the check boxes to select citations from your search results or Clipboard. You may move to other pages and continue your selections. You may also choose to email all citations shown on the page without making any selections, or email an individual item from its abstract page.
2. Click the Email button. - If you are not signed in to My NCBI, you will be redirected to sign in or register for a free My NCBI account if you do not already have one.
- After signing in, you will be redirected back to the email form in PubMed.
3. (Optional) Customize the subject line and enter any text you want to include in the email.
4. Choose a Selection and Format: - Selection: The citations you would like to email. - Selection: The number of selected items will be shown, for example: Selection (87).
- All results on this page
- All results (up to a maximum of 1,000 citations)
- Format: Summary, Summary (text), Abstract, or Abstract (text)
5. Click Send email. The system returns you to your current page and displays a confirmation message.
 More information about emailing citations:
 - Your citations will be sent from the NCBI automatic mail server with the sender's email address [nobody@ncbi.nlm.nih.gov]. Do not reply to this message, as this is not a functioning customer service email address and is not monitored.
- To use the email feature, [register for a free My NCBI account](https://account.ncbi.nlm.nih.gov/signup/). See the [My NCBI Help](https://www.ncbi.nlm.nih.gov/books/NBK3842/) to learn more about My NCBI and the additional features it supports for NCBI databases.
 ### Create an email alert for a search
 Click "Create alert" under the search bar to create an automatic email update for searches. You must sign in to My NCBI to use this feature. See [Saving and Managing Searches](https://www.ncbi.nlm.nih.gov/books/NBK53592/) for more information.
 ### Create an RSS feed for a search
 Click on Create RSS under the search box at the top of the page to create an RSS feed for your search.
 1. The RSS feed name will default to the search terms. You can edit the RSS feed name as needed.
2. Use the pull-down menu to select the number of items displayed. You may manually edit the limit= parameter in the RSS feed link created in Step 4 to display up to a maximum of 200 items. Please note that increasing this limit will also increase the loading time.
3. Click the Create RSS button.
4. The RSS Feed Link will appear; click on Copy to copy the link.
5. Use this link with your feed reader or other application.
 ### Print your search results
 Use the print function of your web browser. To print citations from different searches, save the citations in [PubMed’s Clipboard](#save-citations-to-clipboard), and then print.
 See also:
 - [Showing more results](#show-more-results)
- [Changing the display format](#changing-display-format)
 ### Get a permalink to bookmark or share your search
 To get the URL for an individual citation:
 - From an abstract page, click the Permalink button on the right sidebar and click "Copy".
- From the search results, right click an article title and choose the copy link option; exact wording and functionality may change depending on your operating system and browser.
- From a mobile device, use your device’s native share functionality to copy a link.
 To get the URL for your search results, copy the URL from your web browser's address bar or bookmark the URL using your web browser's bookmark function.
 To create a URL manually:
 1. Use the base URL: https://pubmed.ncbi.nlm.nih.gov/?term=search
2. Replace “search” in the base URL with your query terms
3. Escape spaces by converting them to plus signs (+); for example, Biochem Soc Trans should be entered as: Biochem+Soc+Trans
 The number of characters you can use may be limited by your browser’s maximum URL length (which may be different for each browser).
 Example
 Search PubMed for articles about antioxidant and chocolate:
 [https://pubmed.ncbi.nlm.nih.gov/?term=antioxidant+chocolate](https://pubmed.ncbi.nlm.nih.gov/?term=antioxidant+chocolate)
 Optional search parameters:
 - format=summary, abstract, pubmed, pmid
- sort=relevance, date, pubdate, fauth, jour
- sort_order=asc
- size=10, 20, 50, 100, 200
 Example
 Search PubMed for articles about breast cancer, sorted by ascending publication date (oldest to newest), and display 50 citations per page:
 [https://pubmed.ncbi.nlm.nih.gov/?term=breast+cancer&sort=pubdate&sort_order=asc&size=50](https://pubmed.ncbi.nlm.nih.gov/?term=breast+cancer&sort=pubdate&sort_order=asc&size=50)
 More information about PubMed links:
 - Some settings in PubMed rely on cookies and other session data that may not be present in the URL. For example, searches that were created using a search number in Advanced History (e.g., #1 OR #2 AND human[mh]) cannot be saved using the URL because the search will be lost when your History expires.
- Users intending to send frequent queries or retrieve large numbers of records from the NCBI databases should use [E-Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/). Users must comply with the [usage guidelines and requirements](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen) to prevent overloading NCBI systems.
- The [NCBI Disclaimer and Copyright](https://www.ncbi.nlm.nih.gov/About/disclaimer.html) notice must be evident to users. Users are advised to consult legal counsel to ensure compliance with intellectual property laws. NLM cannot provide advice about copyright issues.
 ### Download PubMed data
 Once a year, NLM releases a complete (baseline) set of PubMed citation records in XML format for download from our FTP servers. Incremental update files are released daily and include new, revised, and deleted citations. The PubMed DTD states any changes to the structure and allowed elements from year to year.
 Note: Binary mode must be used when downloading data from our FTP servers.
 - Documentation: [PubMed XML Elements and Attributes](https://dtd.nlm.nih.gov/ncbi/pubmed/doc/out/250101/index.html)
- [Terms and Conditions](https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/README.txt)
- [PubMed Baseline](https://ftp.ncbi.nlm.nih.gov/pubmed/baseline)
- [PubMed Update Files](https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles)
- [PubMed DTD](https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_250101.dtd)
 For more information, please see [Download PubMed Data](https://pubmed.ncbi.nlm.nih.gov/download/).
  ## Advanced Search
 - [Searching in a specific field](#searching-specific-field)
- [Browsing the index of terms](#browsing-show-index)
- [History](#pubmed-search-history)
- [Previewing the number of search results](#previewing-number-of-search-results)
- [Combining searches using History](#combining-history)
- [Viewing the Search Details](#viewing-search-details)
 Tools included on the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page help users to: search for terms in a specific field, combine searches and build large, complex search strings, see how each query was translated by PubMed, and compare number of results for different queries.
 ### Searching in a specific field
 Use the Advanced Search Builder to search for terms in a specific field, such as author or journal. For some fields, an autocomplete feature will provide suggestions as you type.
 1. From the "All Fields" drop-down menu, select the field you would like to search.
2. Add terms from the builder to the query box to construct your search. The default Boolean operator is AND; if desired, choose OR or NOT from the pull-down menu.
3. Once you have finished adding terms to the query box, click Search (or Add to History) to run the search.
 You may also search a specific field -- and bypass [Automatic Term Mapping](#automatic-term-mapping) -- by adding a [search field tag](#using-search-field-tags) to a term.
 ### Browsing the index of terms
 The Advanced Search Builder includes the Show Index feature, which provides an alphabetical display of terms appearing in selected PubMed search fields. You can browse by all fields or within specific fields such as MeSH Terms.
 1. Click Advanced to navigate to the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page, and use the Builder to select a search field from the All Fields menu. Note: Show Index is not available for every search field. The Show Index link will only display for fields that are compatible with this feature.
2. Enter a term in the search box, then click Show Index.
3. The index displays an alphabetic list of search terms and the approximate number of citations for each term (the actual citation count is returned when the search is executed).
4. Scroll until you find a term you want to include in your search, and then highlight it to add it to the search box.
5. Multiple terms may be selected from the list and added to the search box.
6. Add terms from the builder to the query box to construct your search.
7. Once you have finished adding terms to the query box, click Search (or Add to History) to run the search.
 More information about using the index:
 - PubMed processes all Boolean operators left to right.
- The builder will automatically OR (and add parentheses) for multiple terms selected from the index.
- A slash will display after a space. For example, the MeSH Term and Subheading "zika virus/analysis" will display after "zika virus infection/virology." Enter MeSH terms followed by a slash to go directly to the display for the MeSH/Subheading combination counts in the index.
- Show Index is not available for date fields.
 ### History
 Your PubMed search history appears on the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page under History. This feature requires your web browser to accept cookies.
 Descriptions of each column in the History table appear below:
 - Search: Searches are numbered in chronological order. - Search numbers may be used in place of the search string itself when combining queries (e.g., #1 OR #2).
- A repeated query will move to the top of History but will retain its original numbering.
- History is limited to the last 100 searches. Once the maximum number is reached, PubMed will remove the oldest search from history and add the most current search.
- Actions: Add, delete, or save a query. Adding queries from History places the search string into the Query box to be used in the next search. Deleting a query removes it from History.
- Query: This column shows previous search strings as entered by the user.
- Details: PubMed may modify or add search terms to a search to optimize retrieval, e.g., using automatic term mapping. Click the chevron icon " > " to expand search details and see how the search was translated.
- Results: The total number of citations retrieved for that query. Click the number to run the search and see the results in PubMed.
- Time: Timestamp of when the search was conducted.
- Download: Click Download to generate a CSV file of current History items. - Please note, Microsoft Excel is typically unable to display or print more than a maximum of 1024 characters in a cell; therefore, you may want to open the CSV file with a text editor to display your complete searches.
- Delete: Click "Delete" to remove all queries from History; otherwise, History expires after 8 hours of inactivity.
 ### Previewing the number of search results
 1. Click Advanced to navigate to the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page.
2. Use the builder to add search terms to the query box, or type your search directly into the query box.
3. Use the split button to toggle the button function from "Search" to "Add to History".
4. Click Add to History. This will run the search without leaving the Advanced Search page.
5. See your query including the number of results in the History table.
 ### Combining searches using History
 Searches can be combined or used in later searches using your search History.
 1. Click Advanced to navigate to the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page.
2. In the History table, click the More Actions icon " ... " next to your query.
3. From the available options, select "Add query" to copy the query to the Query box.
4. After you've added content to the Query box, options to use the Boolean operators AND, OR, or NOT will appear when adding more queries to the Query box.
5. Edit your query in the Query box if you would like to make any changes before running the search.
6. Click Search (or Add to History).
 More information about combining searches from your History:
 - Search numbers may be used in place of the search string itself when combining queries (e.g., #1 OR #2).
- Citations in the Clipboard are represented by the search number #0, which may be used in searches. For example, to limit the citations you have collected in the clipboard to English language citations, use the following search: #0 AND english [la]. This does not change or replace the Clipboard contents.
 ### Viewing the Search Details
 PubMed may modify or add additional search terms to your search to optimize retrieval, such as: MeSH terms, British/American spellings, singular/plural word forms, and other synonyms.
 1. Search Details are included on the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page under History.
2. Click the chevron icon " > " next to a query in History to expand the Search Details.
3. When expanded, the details below a query in the History table show the search strategy used to run the search.
 More information about search details:
 - Translations show individual term mappings using PubMed's search rules and syntax. Query terms without translations will not be listed in this section; for example, [exact phrases](#searching-for-a-phrase) bypass [Automatic Term Mapping (ATM)](#automatic-term-mapping).
- Warnings are displayed for the original query with potential errors in bold and red type, such as syntax errors, terms not found, or invalid tags. Warnings also appear as a highlighted message in PubMed when the search is run or added to History.
  ## Other services
 - [Clinical Queries](#clinical-queries)
- [Single Citation Matcher](#single-citation-matcher)
- [Search PubMed using the MeSH database](#using-mesh-database)
- [Search for journal information in the NLM Catalog](#searching-nlm-catalog)
- [Using the E-utilities API tools](#e-utils-api)
- [Citation Matcher API](#citation-matcher-api)
- [Batch Citation Matcher](#batch-citation-matcher)
- [Consumer health](#consumer-health)
 ### Clinical Queries
 [PubMed Clinical Queries](https://pubmed.ncbi.nlm.nih.gov/clinical/) provides specialized searches for:
 - [COVID-19 Articles](#covid19-articles)
- [Clinical Study Categories](#clinical-study-category)
- [Medical Genetics](#medical-genetics)
 #### Search for COVID-19 articles
 The COVID-19 article filters limit retrieval to citations about the 2019 novel coronavirus. Results are displayed in a column filtered by research topic categories. See [COVID-19 article filters](#covid19-article-filters) for the filter search strategies; these may evolve over time.
 To find citations using the COVID-19 article filters:
 1. Click Clinical Queries from the PubMed homepage
2. Enter your search terms in the search box
3. Click Search
4. Select a Category: General, Mechanism, Transmission, Diagnosis, Treatment, Prevention, Case Report, Forecasting, or Long COVID
5. Preview results in the COVID-19 Articles column
6. To view the results in PubMed, click the "See all" link below the results preview
 To use the COVID-19 article filters in a query, add the filter name to your search with the search field tag [Filter], e.g., LitCPrevention[Filter]. The available filters are:
 - LitCGeneral
- LitCMechanism
- LitCTransmission
- LitCDiagnosis
- LitCTreatment
- LitCPrevention
- LitCCaseReport
- LitCForecasting
- LitCLongCOVID
 Example
 Search PubMed for Remdesivir with the COVID-19 General filter:
 [Remdesivir AND LitCGeneral[Filter]](https://pubmed.ncbi.nlm.nih.gov/?term=Remdesivir+AND+LitCGeneral%5Bsb%5D)
 #### Search by clinical study category
 Clinical Study Categories use a specialized search method with built-in search filters that limit retrieval to citations reporting research conducted with specific methodologies, including those that report applied clinical research. See [Clinical Study Categories filters](#clinical-study-category-filters) for the filter search strategies.
 To find citations using the Clinical Study Categories:
 1. Click Clinical Queries from the PubMed homepage
2. Enter your search terms in the search box
3. Click Search
4. Select a Category: Therapy, Diagnosis, Etiology, Prognosis, or Clinical Prediction Guides
5. Select a Scope: Narrow (specific search) or Broad (sensitive search)
6. Preview results in the Clinical Study Categories column
7. To view the results in PubMed, click the "See all" link below the results preview
 #### Medical genetics searches
 The Medical Genetics filters limit retrieval to citations related to various topics in medical genetics. See [Medical genetics search filters](#medical-genetics-filters) for the filter search strategies.
 To use a Medical Genetics filter, add the filter name to your search with the search field tag [Filter], e.g., Genetic Testing[Filter]. The available filters are:
 - Diagnosis
- Differential Diagnosis
- Clinical Description
- Management
- Genetic Counseling
- Molecular Genetics
- Genetic Testing
- Medical Genetics
 Example
 Search PubMed for sickle cell anemia using the Genetic Counseling filter:
 [sickle cell anemia AND genetic counseling[Filter]](https://pubmed.ncbi.nlm.nih.gov/?term=sickle+cell+anemia+AND+genetic+counseling%5Bfilter%5D)
 ### Single Citation Matcher
 The [Single Citation Matcher](https://pubmed.ncbi.nlm.nih.gov/citmatch/) has a fill-in-the-blank form for searching for a citation when you have some bibliographic information, such as journal name, volume, or page number.
 - Click Single Citation Matcher on the PubMed homepage.
- Enter the citation information.
- Click Go.
 More information about using the Single Citation Matcher:
 - The journal box includes an autocomplete feature that suggests titles as you enter a title abbreviation or full title. Titles displayed by the autocomplete menu are in ranked order based on the number of citations in PubMed.
- After selecting a journal with special characters (e.g., ampersand, colon) when using the Back button to return to the Single Citation Matcher you must clear and reenter the title.
- The author box also includes an autocomplete feature that suggests author names in ranked order based on the number of citations. Full author names may be searched for citations published from 2002 forward if the full author name is available in the article.
- Click either the 'Only as first author' or ‘Only as last author’ check box to limit an author name to the first or last author.
 ### Search PubMed using the MeSH database
 [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) (Medical Subject Headings) is the NLM controlled vocabulary thesaurus used for indexing PubMed citations.
 Use the [MeSH database](https://www.ncbi.nlm.nih.gov/mesh/) to find MeSH terms, including Subheadings, Publication Types, Supplementary Concepts and Pharmacological Actions, and then build a PubMed search. The MeSH database can be searched by MeSH term, MeSH Entry Term, Subheading, Publication Type, Supplementary Concept, or MeSH Scope Note.
 More information about the MeSH database:
 - An autocomplete feature is available from the search box.
- Search results are displayed in relevance-ranked order, therefore, when a user’s search exactly matches a MeSH Term, that Term is displayed first.
- Click the MeSH term from the Summary display or choose Full from the display format menu to view additional information and search specifications, such as Subheadings, restrict to Major MeSH Topic, or exclude terms below the term in the MeSH hierarchy.
- Year Introduced is the year the term was added to MeSH. If more than one year is shown, the term was available for indexing back to the earliest year noted. Articles are indexed using the vocabulary in place at the time of indexing, therefore, the year introduced for a term and the date of publication of a citation indexed with that term may not agree.
 #### Launch PubMed searches from the MeSH database
 To build a PubMed search from MeSH:
 1. Run a search in the [MeSH database](https://www.ncbi.nlm.nih.gov/mesh/).
2. Select terms using the check boxes.
3. Click "Add to search builder" in the PubMed search builder portlet.
4. You may continue searching and including additional terms to the PubMed search builder using the "Add to search builder" and Boolean pull-down menu.
5. When you are finished, click "Search PubMed."
 ### Search for journal information in the NLM Catalog
 The [NLM Catalog](https://www.ncbi.nlm.nih.gov/nlmcatalog) includes information about the journals in PubMed and the other NCBI databases.
 Click [Journals in NCBI Databases](https://www.ncbi.nlm.nih.gov/nlmcatalog/journals) on the homepage of NLM Catalog or the [Journals](https://www.ncbi.nlm.nih.gov/nlmcatalog/journals) link on the PubMed homepage to limit your NLM Catalog results to the subset of journals that are referenced in NCBI database records.
 See the [NLM Catalog help](https://www.ncbi.nlm.nih.gov/books/n/helpcatalog/catalog/#catalog.Searching_for_Journals_in_NLM__1) for additional information.
 Other journal resources include:
 - PubMed journals with [links to full text](https://www.ncbi.nlm.nih.gov/projects/linkout/journals/jourlists.fcgi)
- List of [all journals included in PubMed via FTP](#journal-lists)
- [List of Serials](https://www.nlm.nih.gov/tsd/serials/lsiou.html) Indexed for Online Users
 ### Using the E-utilities API tools
 [E-utilities](https://www.ncbi.nlm.nih.gov/books/n/helpeutils/chapter1/) are tools that provide access to data outside of the regular NCBI web search interface. This may be helpful for retrieving search results for use in another environment. If you are interested in large-scale data mining on PubMed data, you may download the data for free from our [FTP server](https://ftp.ncbi.nlm.nih.gov/pubmed/baseline). Please see the [terms and conditions](https://www.nlm.nih.gov/databases/download/terms_and_conditions.html) for data users.
 ### Citation Matcher API
 - [API root](#citation-matcher-api-root)
- [Fielded search](#citation-matcher-fielded-search)
- [Heuristic search](#citation-matcher-heuristic-search)
- [Auto search](#citation-matcher-auto-search)
- [Rate control](#citation-matcher-rate-control)
 The PubMed Citation Matcher API finds PubMed identifiers (PMIDs) for citation data in structured or raw form. The interface supports three retrieval methods:
 1. `field` - runs a fielded search using core bibliographic information, such as journal, date, or volume.
2. `heuristic` - collects all input elements into a single string and returns the closest matching documents.
3. `auto` - combines the two above methods and switches to heuristic mode if the fielded search has not yielded a result. This is the default method.
 More information about the Citation Matcher API:
 - The API supports both GET and POST requests.
- Data is exchanged in JSON.
- Input data should be UTF-8 encoded.
- The API returns a maximum of 20 PMIDs; queries returning more than 20 PMIDs are treated as bad requests.
 #### API root
 The API root is:
 ```
https://pubmed.ncbi.nlm.nih.gov/api/citmatch/
```
 #### Fielded search
 `method=field` runs a fielded search using core bibliographic information, such as journal, date, or volume. This functionality is similar to [E-utilities ESearch](https://www.ncbi.nlm.nih.gov/books/n/helpeutils/chapter1/); users should select the API that best suits their needs.
 For a structured search, the following fields can be used:
 - `journal` - the name of the journal
- `pdat` - the publication date, in the format YYYY/MM/DD
- `volume` - the volume of the publication
- `issue` - the volume of the publication
- `authors` - one or more author names, in the format "Surname Initial" (Doe J). Optionally, the position may be specified as first, last, or auto.
 Example
 Example fielded search:
 [https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=field&journal=Front%20Immunol&volume=13&page=826091](https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=field&journal=Front%20Immunol&volume=13&page=826091)
 Example
 GET request URL:
 ```
/citmatch/?method=field&journal=Front+Immunol&volume=13&page=826091&authors=Ihunwo+A
```
 POST request data:
 ```
{
    "citmatch": {
        "method": "field",
        "journal": "Front Immunol",
        "volume": "13",
        "page": "826091",
        "authors": [
            {
                "name": "Ihunwo AO",
                "position": "first"
            }
        ]
    }
}
            
```
 Response:
 ```
{
    "version": "1.0",
    "operation": "citmatch",
    "success": true,
    "result": {
        "count": 1,
        "type": "uids",
        "uids": [
            {"pubmed": "35251006"}
        ]
    }
}                
            
```
 #### Heuristic search
 `method=heuristic` collects all input elements into a single string and returns the closest matching documents. It is sufficient to supply a raw citation string, such as: "The role of drag in insect hovering. J. Exp. Biol. 2004;207:4147–4155."
 Example
 Example heuristic search:
 [https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=heuristic&raw-text=Neurogenesis%20and%20Viral%20Infection.%20Front%20Immunol.%202022%20Feb%2017;13:82609](https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=heuristic&raw-text=Neurogenesis%20and%20Viral%20Infection.%20Front%20Immunol.%202022%20Feb%2017;13:82609)
 Example
 GET request URL:
 ```
/citmatch/?method=heuristic&raw-text=Neurogenesis+and+Viral+Infection.+Front+Immunol.+2022+Feb+17;13:82609.
```
 POST request data:
 ```
{
    "citmatch": {
        "method": "heuristic",
        "raw-text": {"text": "Neurogenesis+and+Viral+Infection.+Front+Immunol.+2022+Feb+17;13:82609."}
    }
}
            
```
 Response:
 ```
{
    "version": "1.0",
    "operation": "citmatch",
    "success": true,
    "result": {
        "count": 1,
        "type": "uids",
        "uids": [
            {"pubmed": "35251006"}
        ]
    }
}              
            
```
 #### Auto search
 `method=auto` first runs a [fielded search](#citation-matcher-fielded-search), and if no results are found, it combines the fields and runs a [heuristic search](#citation-matcher-heuristic-search). This is the default method.
 For a structured search, the following fields can be used:
 - `journal` - the name of the journal
- `pdat` - the publication date, in the format YYYY/MM/DD
- `volume` - the volume of the publication
- `issue` - the volume of the publication
- `authors` - one or more author names, in the format "Surname Initial" (Doe J). Optionally, the position may be specified as first, last, or auto.
 Example
 Example auto search:
 [https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=auto&journal=Front+Immunol&volume=13&page=826091&authors=Ranganathan+S](https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=auto&journal=Front+Immunol&volume=13&page=826091&authors=Ranganathan+S) First, a fielded search is run but no results are found due to the specified author not appearing on the citation:
 [https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=field&journal=Front+Immunol&volume=13&page=826091&authors=Ranganathan+S](https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=field&journal=Front+Immunol&volume=13&page=826091&authors=Ranganathan+S) Then it concatenates the fields and runs a heuristic search, which returns the closest matching document(s):
 [https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=heuristic&raw-text=Front+Immunol+13+826091+Ranganathan+S](https://pubmed.ncbi.nlm.nih.gov/api/citmatch/?method=heuristic&raw-text=Front+Immunol+13+826091+Ranganathan+S)
 #### Rate control
 When using the PubMed Citation Matcher API programmatically, we request that you limit your application's rate to 3 requests / sec and do not make concurrent requests to this service, even at off-peak times. Additionally, requests must contain the name of the calling project in the User-Agent HTTP header value; e.g. `Hydra/1.3.15`.
 ### Batch Citation Matcher
 Use the [Batch Citation Matcher](https://pubmed.ncbi.nlm.nih.gov/batchcitmatch/) to retrieve PMIDs for multiple citations. The Batch Citation Matcher requires that you enter the bibliographic information (journal, volume, page, etc.) in a specific format.
 To retrieve PubMed PMIDs:
 1. Create citation strings for the items you would like to retrieve using the following format: 

journal_title|year|volume|first_page|author_name|your_key|


 Fields must be separated by a vertical bar with a final bar at the end of the string.
2. Enter your email address. Email messages may take several minutes to process and be sent to your email address.
3. Upload your citation strings as a text file (.txt) or enter each citation string on a separate line in the text box. If citation strings are entered in the text box and a file is uploaded, the results will be an aggregate of both.
4. Click search.
 If a match is not found the citation string will display one of the following:
 - your_key|NOT_FOUND;INVALID_JOURNAL - The journal name is not valid. See the [journal lists](#journal-lists) or the [NLM Catalog](https://www.ncbi.nlm.nih.gov/nlmcatalog/journals/) to find the correct journal abbreviation.
- NOT_FOUND - The journal name is valid, but the citation string did not find a match.
- AMBIGUOUS - The information provided matches more than one citation. Citation information with 3 or fewer matches include the PMIDs, and more than 3 matches include the total PMID match count. Use the [Single Citation Matcher](https://pubmed.ncbi.nlm.nih.gov/citmatch/) or [ESearch](https://www.ncbi.nlm.nih.gov/books/n/helpeutils/chapter1/) to retrieve all citations for searched fields.
 Notes:
 - Text (.txt) format must be used when uploading a file.
- You may receive multiple emails for searches containing more than 2,000 citation strings.
- Enter author names without punctuation as smith jc. Initials are optional.
- Your key is any string you choose to tag the citation, it is returned unaltered.
- The journal title field may include the full journal title or the NLM title abbreviation.
- Each citation field is searched starting with the journal title until a unique match is found.
- The journal title is a required field however you may omit other fields. If you omit fields you must retain the vertical bars in the citation string. For example, if you omit the volume number 88 from the first example below it should be entered as:


proc natl acad sci u s a|1991||3248|mann bj|P32022-1|
 Example input:
 - proc natl acad sci u s a|1991|88|3248|mann bj|P32022-1|
- proc natl acad sci u s a|1992|89|3271|gould se|P26261-1|
- proc natl acad sci u s a|1970|89|3271|smith|P26261-1|
- res microbiol|1992|143|467|ivey dm|P25966-1|
- science|1987|235|182|palmenberg ac|P12296-2|
- eschatology|1993|12|22|public jq|C12233-2|
- virology|1993|193|492|hardy me|Q02945-1|
- virus genes|1992|6|393||P27423-1|
- yeast|1992|8|253|sasnauskas k|P24813-1|
 Example output:
 - proc natl acad sci u s a|1991|88|3248|mann bj|P32022-1|2014248
- proc natl acad sci u s a|1992|89|3271|gould se|P26261-1|1565618
- proc natl acad sci u s a|1970|89|3271|smith|P26261-1|NOT_FOUND
- res microbiol|1992|143|467|ivey dm|P25966-1|1448623
- science|1987|235|182|palmenberg ac|P12296-2|3026048
- C12233-2|NOT_FOUND;INVALID_JOURNAL
- virology|1993|193|492|hardy me|Q02945-1|8382410
- virus genes|1992|6|393||P27423-1|1335631
- yeast|1992|8|253|sasnauskas k|P24813-1|1514324
 ### Consumer health
 The National Library of Medicine cannot provide specific medical advice. NLM urges you to consult a qualified health care professional for answers to your medical questions. NLM does not have pamphlets or other materials to mail.
 [MedlinePlus](https://www.nlm.nih.gov/medlineplus/) and MedlinePlus en español are specifically designed for consumers, containing hundreds of topic pages including NIH-written descriptive information, videos, health check tools, drug, herb and supplement info, links to Fact Sheets from other NIH Institutes, the CDC, etc., and more.
  ## Appendices
 - [Further assistance and training](#assistance-training)
- [How PubMed works: Automatic Term Mapping (ATM)](#automatic-term-mapping)
- [Algorithm for finding best matching citations in PubMed](#best-match-algorithm)
- [PubMed coverage](#pubmed-coverage)
- [PubMed format](#pubmed-format)
- [PubMed data field descriptions](#data-fields)
- [NLM author indexing policy](#author-indexing-policy)
- [Error messages](#error-messages)
- [Cookies](#help-cookies)
- [MeSH Subheadings](#mesh-subheadings)
- [Stopwords](#help-stopwords)
- [PubMed character conversions](#character-conversions)
- [Publication Types](#publication-types)
- [Status Subsets](#citation-status-subsets)
- [Filter search strategies](#filter-strategies)
- [Clinical Queries filters](#clinical-queries-filters)
- [Computation of similar articles](#computation-of-similar-articles)
- [Journal lists](#journal-lists)
 ### Further assistance and training
 #### Contact customer support
 - E-mail the [PubMed Help Desk](https://support.nlm.nih.gov/support/create-case/?category=pubmed)
- Call the NLM Customer service desk: 1-888-FIND-NLM (1-888-346-3656)
 #### Other NLM publications
 - [PubMed Online Training](https://learn.nlm.nih.gov/documentation/training-packets/T0042010P/)
- [PubMed Trainer's Toolkit](https://learn.nlm.nih.gov/documentation/training-packets/T0022014P/)
- [NLM Technical Bulletin](https://www.nlm.nih.gov/pubs/techbull/tb.html)
 ### How PubMed works: Automatic Term Mapping (ATM)
 Untagged terms that are entered in the search box are matched (in this order) against a Subject translation table (including [MeSH (Medical Subject Headings)](https://www.nlm.nih.gov/mesh/meshhome.html)), a Journals translation table, the Author index, and an Investigator (Collaborator) index.
 When a match is found for a term or phrase in a translation table the mapping process is complete and does not continue on to the next translation table.
 To see how your terms were translated, check the [Search Details](#viewing-search-details) available on the [Advanced Search](https://pubmed.ncbi.nlm.nih.gov/advanced/) page for each query under History. If you want to report a translation that does not seem accurate for your search topic, please e-mail the information to the [NLM Help Desk](https://support.nlm.nih.gov/support/create-case/?category=pubmed).
 #### 1. Subject translation table
 The Subject Translation Table contains:
 - British and American spellings
- Pairs: singular and plural word forms, synonyms, and other closely related terms
- Drug brand name to generic name translations
- [MeSH](https://www.nlm.nih.gov/mesh/meshhome.html) terms
- The See-Reference mappings (also known as entry terms) for MeSH terms
- [MeSH Subheadings](#mesh-subheadings)
- [Publication Types](#publication-types)
- [Pharmacologic action](https://www.ncbi.nlm.nih.gov/mesh/68020228) terms
- Terms derived from the [Unified Medical Language System (UMLS)](https://www.nlm.nih.gov/research/umls/) that have equivalent synonyms or lexical variants in English
- Supplementary concept (substance) names and their synonyms.
 If a match is found in this translation table, the term will be searched as MeSH (that includes the MeSH term and any specific terms indented under that term in the MeSH hierarchy), and in all fields.
 For example, if you enter child rearing in the search box, PubMed will translate this search to: "child rearing"[MeSH Terms] OR ("child"[All Fields] AND "rearing"[All Fields]) OR "child rearing"[All Fields]
 If you enter a MeSH Term that is also a Pharmacologic Action PubMed will search the term as [MeSH Terms], [Pharmacologic Action], and [All Fields].
 If you enter an entry term for a MeSH term the translation will also include an all fields search for the MeSH term associated with the entry term. For example, a search for odontalgia will translate to: "toothache"[MeSH Terms] OR "toothache"[All Fields] OR "odontalgia"[All Fields] OR "odontalgias"[All Fields] because Odontalgia is an entry term for the MeSH term toothache.
 Substance name mappings do not include a mapping for individual terms in a phrase, e.g., IL-22 will not include IL[All Fields] AND 22[All Fields].
 MeSH term mappings that include a standalone number or single character do not include a mapping for individual terms in a phrase, e.g., Protein C will not include Protein[All Fields] or C[All Fields].
 #### 2. Journals translation table
 The Journals translation table contains the:
 - full journal title
- title abbreviation
- ISSN and eISSN number.
 These will automatically map to the journal abbreviation that is used to search journals in PubMed and in all fields. For example, a search for endocrine pathology will translate to: "Endocr Pathol"[Journal] OR ("endocrine"[All Fields] AND "pathology"[All Fields]) OR "endocrine pathology"[All Fields]
 #### 3. Author index
 If the term is not found in the above tables, and is not a single term, PubMed checks the author index for a match. The author index includes author names and initials, as well as full author names for articles published from 2002 forward, if available.
 More information about author searching:
 - PubMed automatically truncates a search for an author's name to account for varying initials, e.g., o'brien j retrieves o'brien ja, o'brien jb, o'brien jc jr, as well as o'brien j.
- When combining multiple authors, to avoid a match with full author names, include initials or use the [au] search tag, e.g., ryan[au] james[au]. Author names comprised of only stopwords, e.g., as a, are not searched as authors if they are part of phrase, chemical burn as a danger, unless the search only includes the author name, e.g., as a.
- Enter a full author name in natural or inverted order, e.g., julia s wong or wong julia s.
- A comma following the last name for searching is optional. For some names, however, it is necessary to distinguish which name is the last name by using the comma following the last name, e.g., james, ryan.
- Omit periods after initials and put all suffixes at the end, e.g., vollmer charles jr
- Initials and suffixes are not required, if you include a middle initial or suffix, you will only retrieve citations for articles that were published using the middle initial or suffix.
- To distinguish author initials that may match a full author name use the [fau] search tag, e.g., peterson do[fau].
 #### 4. Investigator (Collaborator) index
 If the term is not found in the above tables, except for Author, and is not a single term, the investigator index is consulted for a match. The investigator (collaborator) index includes full names, if available. Enter a full investigator name in natural or inverted order, e.g., harry janes or janes harry.
 #### 5. If no match is found?
 PubMed breaks apart the phrase and repeats the above automatic term mapping process until a match is found. PubMed ignores [stopwords](#help-stopwords) in searches.
 If there is no match, the individual terms will be combined (ANDed) together and searched in all fields.
 When a search includes terms that were tagged with a search field during the automatic term mapping process and retrieves zero results, the system triggers a subsequent search using "Schema: all ." "Schema: all" modifies the search by removing the automatically added search field tags, and then searches each term in all fields.
 ### Algorithm for finding best matching citations in PubMed
 The learned ranking algorithm combines over 150 signals that are helpful for finding best matching results. Most of these signals are computed from the query-document term pairs (e.g., number of term matches between the query and the document) while others are either specific to a document (e.g., publication type; publication year) or query (e.g., query length). The new ranking model was built on relevance data extracted from the anonymous and aggregated PubMed search logs over an extended period of time.
 For more information about the Best Match algorithm, please see:
 - Technical details in the paper [Best Match: New relevance search for PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6112631/) by Fiorini N, Canese K, Starchenko G, et al. in PLoS Biol (2018).
- NLM Technical Bulletin article: [Updated Algorithm for the PubMed Best Match Sort Order](https://www.nlm.nih.gov/pubs/techbull/jf17/jf17_pm_best_match_sort.html)
 ### PubMed coverage
 The PubMed database contains citations and abstracts to biomedical literature, facilitating searching across several NLM literature resources:
 - [MEDLINE](#coverage-medline)
- [PubMed Central (PMC)](#coverage-pmc)
- [NCBI Bookshelf](#coverage-bookshelf)
 For additional information, please see the webpages [About PubMed](https://pubmed.ncbi.nlm.nih.gov/about/) and [About NLM's Literature Databases](https://www.nlm.nih.gov/portals/publishers.html).
 PubMed includes citations to original research articles, literature reviews, case reports, letters, editorials, commentaries, and other selected publications on scientific and medical topics (see: [publication types found in PubMed](#publication-types)). Some categories of content are out of scope for PubMed, such as: book reviews, individual conference abstracts, [obituaries and in memoriam articles](https://www.nlm.nih.gov/pubs/techbull/mj18/brief/mj18_obituaries_in_memoriam.html), news and announcements, and brief summaries of research articles. More examples are included in [XML Help for PubMed Data Providers: What types of articles are accepted?](https://www.ncbi.nlm.nih.gov/books/NBK3828/#publisherhelp.What_types_of_articles_are).
 #### MEDLINE
 [MEDLINE](https://www.nlm.nih.gov/bsd/medline.html) contains citations to journal articles in the life sciences with a concentration on biomedicine. The MEDLINE database contains citations from the late [1940s to the present](https://www.nlm.nih.gov/databases/databases_oldmedline.html), with some older material.
 New citations from MEDLINE journals are received electronically from publishers and appear in PubMed daily. Most citations progress to in-process, and then to indexed for MEDLINE; however, not all citations will be indexed for MEDLINE. PubMed includes some citations from MEDLINE journals that are not indexed for MEDLINE, such as:
 - Citations preceding the date that a journal was selected for MEDLINE indexing.
- Out-of-scope citations (e.g., articles on plate tectonics or astrophysics) from certain MEDLINE journals, primarily general science and chemistry journals, for which the life sciences articles are indexed for MEDLINE.
 Citations that have been indexed for MEDLINE and updated with NLM [Medical Subject Headings (MeSH)](https://www.nlm.nih.gov/mesh/meshhome.html), publication types, GenBank accession numbers, and other indexing data are available daily. To limit your search to MEDLINE citations, add medline[sb] to your search.
 ##### Indexing method
 The method used to assign Medical Subject Headings (MeSH) has changed over time. For more information, please see [Incorporating Values for Indexing Method in MEDLINE/PubMed XML](https://www.nlm.nih.gov/pubs/techbull/ja18/ja18_indexing_method.html). Use the following searches to find citations indexed with each method:
 - Automated - MeSH indexing is provided algorithmically. Search: [indexingmethod_automated](https://pubmed.ncbi.nlm.nih.gov/?term=indexingmethod_automated)
- Curated - MeSH indexing is provided algorithmically and a human reviewed (and possibly modified) the algorithm results. Search: [indexingmethod_curated](https://pubmed.ncbi.nlm.nih.gov/?term=indexingmethod_curated)
- Manual – MeSH indexing is provided by humans with no algorithmic suggestions. Search: [indexingmethod_manual](https://pubmed.ncbi.nlm.nih.gov/?term=indexingmethod_manual)
 #### PubMed Central (PMC)
 [PubMed Central (PMC)](https://www.ncbi.nlm.nih.gov/pmc/about/intro/) is a full text archive that includes articles from journals reviewed and selected by NLM for archiving (current and historical), as well as individual articles and [preprints](#coverage-preprints) collected for archiving in compliance with funder policies. Some PMC content is not cited in PubMed, such as book reviews and conference abstracts (see: [PubMed coverage](#pubmed-coverage)).
 ##### Preprints
 As of June 2020, PubMed Central (PMC) includes preprints that report NIH-funded research results. Citations to these preprints are deposited in PubMed. To learn more, see: [NIH Preprint Pilot](https://www.ncbi.nlm.nih.gov/pmc/about/nihpreprints/).
 Example
 To search for preprints in PubMed, include preprint[filter] in your query.
 [covid-19 AND preprint[pt]](https://pubmed.ncbi.nlm.nih.gov/?term=covid-19+AND+preprint%5Bpt%5D)
 Example
 To exclude preprints from your search results in PubMed, use the Boolean operator NOT.
 [covid-19 NOT preprint[pt]](https://pubmed.ncbi.nlm.nih.gov/?term=covid-19+NOT+preprint%5Bpt%5D)
 #### NCBI Bookshelf
 [Bookshelf](https://www.ncbi.nlm.nih.gov/books/) is a full text archive of books, reports, databases, and other documents related to biomedical, health, and life sciences. PubMed includes citations for books and some individual chapters available on Bookshelf.
 ### PubMed format
 The PubMed Format tags table defines the data tags that compose the PubMed format. The tags are presented in alphabetical order. Some of the tags (e.g., CIN) are not mandatory and therefore will not be found in every record. Other tags (e.g., AU, MH, and RN) may occur multiple times in one record. You can download records in PubMed format [as a text file (.txt)](#save-citations-to-file) or as an .nbib file for [exporting into citation management software programs](#citation-management-software).
 Not all fields are searchable in PubMed. See [Search field tags](#search-tags) for the list of searchable fields.
 PubMed Format tags Tag Name Description AB Abstract English language abstract taken directly from the published article AD Affiliation Author or corporate author addresses AID Article Identifier Article ID values supplied by the publisher may include the pii (controlled publisher identifier), doi (digital object identifier), or book accession AU Author Authors AUID Author Identifier Unique identifier associated with an author, corporate author, or investigator name BTI Book Title Book Title CI Copyright Information Copyright statement provided by the publisher CIN Comment In Reference containing a comment about the article CN Corporate Author Corporate author or group names with authorship responsibility COI Conflict of Interest Conflict of interest statement CON Comment On Reference upon which the article comments CP Chapter Book chapter CRDT Create Date The date the citation record was first created CRF Corrected and republished from Final, correct version of an article CRI Corrected and republished in Original article that was republished in corrected form CTDT Contribution Date Book contribution date CTI Collection Title Collection Title DCOM Completion Date NLM internal processing completion date DDIN Dataset described in Citation for the primary article resulting from a dataset DRIN Dataset use reported in Citation for an article that uses a dataset from another scientific article DEP Date of Electronic Publication Electronic publication date DP Publication Date The date the article was published DRDT Date Revised Book Revision Date ECF Expression of Concern For Reference containing an expression of concern for an article ECI Expression of Concern In Cites the original article for which there is an expression of concern EDAT Entry Date The date the citation was added to PubMed; the date is set to the publication date if added more than 1 year after the date published EFR Erratum For Cites the original article for which there is a published erratum; as of 2016, partial retractions are considered errata EIN Erratum In Cites a published erratum to the article ED Editor Book editors EN Edition Book edition FAU Full Author Name Full author names FED Full Editor Name Full editor names FIR Full Investigator Name Full investigator or collaborator names FPS Full Personal Name as Subject Full Personal Name of the subject of the article GN General Note Supplemental or descriptive information related to the document GR Grants and Funding Grant numbers, contract numbers, and intramural research identifiers associated with a publication GS Gene Symbol Abbreviated gene names (used 1991 through 1996) IP Issue The number of the issue, part, or supplement of the journal in which the article was published IR Investigator Investigator or collaborator IRAD Investigator Affiliation Investigator or collaborator addresses IS ISSN International Standard Serial Number of the journal ISBN ISBN International Standard Book Number JID NLM Unique ID Unique journal ID in the NLM catalog of books, journals, and audiovisuals JT Full Journal Title Full journal title from NLM cataloging data LA Language The language in which the article was published LID Location ID The pii or doi that serves the role of pagination LR Modification Date Citation last revision date MH MeSH Terms NLM Medical Subject Headings (MeSH) controlled vocabulary MHDA MeSH Date The date MeSH terms were added to the citation. The MeSH date is the same as the Entrez date until MeSH are added MID Manuscript Identifier Identifier assigned to an author manuscript submitted to the NIH Manuscript Submission System NM Substance Name Supplementary Concept Record (SCR) data OAB Other Abstract Abstract supplied by an NLM collaborating organization OABL Other Abstract Language Language of an abstract available from the publisher OCI Other Copyright Information Copyright owner OID Other ID Identification numbers provided by organizations supplying citation data ORI Original Report In Cites the original article associated with the patient summary OT Other Term Non-MeSH subject terms (keywords) either assigned by an organization identified by the Other Term Owner, or generated by the author and submitted by the publisher OTO Other Term Owner Organization that may have provided the Other Term data OWN Owner Organization acronym that supplied citation data PB Publisher Publishers of Books & Documents citations PG Pagination The full pagination of the article PHST Publication History Status Date Publisher supplied dates regarding the article publishing process and PubMed date stamps: - received: manuscript received for review
- revised: manuscript revised by publisher or author
- accepted: manuscript accepted for publication
- aheadofprint: published electronically prior to final publication
- entrez: PubMed Create Date [crdt]
- pubmed: PubMed Entry Date [edat]
- medline: PubMed MeSH Date [mhda] PL Place of Publication Journal's (country only) or book’s place of publication PMC PubMed Central Identifier Unique identifier for the cited article in PubMed Central (PMC) PMCR PMC Release Availability of PMC article PMID PubMed Unique Identifier Unique number assigned to each PubMed citation PS Personal Name as Subject Individual is the subject of the article PST Publication Status Publication status PT Publication Type The type of material the article represents RF Number of References Number of bibliographic references for Review articles RIN Retraction In Retraction of the article RN EC/RN Number Includes chemical, protocol or disease terms. May also include a number assigned by the Enzyme Commission or by the Chemical Abstracts Service. ROF Retraction Of Article being retracted RPF Republished From Article being cited has been republished or reprinted in either full or abridged form from another source RPI Republished In Article being cited also appears in another source in either full or abridged form RRI Retracted and Republished In Final, republished version of an article RRF Retracted and Republished From Original article that was retracted and republished SB Subset Journal or citation subset values representing specialized topics SFM Space Flight Mission NASA-supplied data space flight/mission name and/or number SI Secondary Source ID Identifies secondary source databanks and accession numbers of molecular sequences discussed in articles SO Source Composite field containing bibliographic information SPIN Summary For Patients In Cites a patient summary article STAT Status Tag Used for internal processing at NLM TA Journal Title Abbreviation Standard journal title abbreviation TI Title The title of the article TT Transliterated Title Title of the article originally published in a non-English language, in that language UIN Update In Update to the article UOF Update Of The article being updated VI Volume Volume number of the journal VTI Volume Title Book Volume Title ### PubMed data field descriptions
 This documentation describes the fields found in PubMed records. If a field is searchable, the search tag appears after the field name in square brackets: Affiliation [ad]. A small number of searchable fields do not correspond to a specific field in the [PubMed format](#pubmed-format).
 - See [Search field tags](#search-tags) for a list of searchable fields.
- See [PubMed format](#pubmed-format) for a quick table view of the fields found in PubMed records.
 #### Affiliation [ad]
 Affiliation may be included for authors, corporate authors and investigators, e.g., cleveland [ad] AND clinic [ad], if submitted by the publisher.
 Multiple affiliations were added to citations starting from 2014, previously only the first author’s affiliation was included. PubMed includes the note "Contributed equally" in the affiliation field when this information is supplied by publishers.
 Searching for terms in the affiliation field searches in all author affiliations on a citation. For example, a search for Hopkins[ad] AND Bloomberg[ad] can find these terms spread across multiple authors’ affiliations on the same citation.
 To search for multiple terms appearing within the same affiliation, use a [proximity search](#proximity-searching). You can also search affiliations using a [phrase search](#searching-for-a-phrase); however, we suggest using a proximity search for more comprehensive results because affiliation data may be provided in a variety of ways for the same institution.
 Example
 Use [proximity searching](#proximity-searching) to find citations with authors from the Johns Hopkins Bloomberg School of Public Health:
 ["Hopkins Bloomberg Public"[ad:~45]](https://pubmed.ncbi.nlm.nih.gov/?term=%22Hopkins+Bloomberg+Public%22%5Bad%3A%7E45%5D) This search will find any citation where the words "Hopkins," "Bloomberg," and "Public" appear in the same affiliation, with no more than forty-five words between each term. Search results may include:
 - Johns Hopkins Bloomberg School of Public Health
- Johns Hopkins University, Bloomberg School of Public Health
- Bloomberg School of Public Health, Johns Hopkins University
- Bloomberg Johns Hopkins University School of Public Health
- ...and more!
 #### All Fields [all]
 Untagged terms and terms tagged with [all] are processed using [Automatic Term Mapping (ATM)](#automatic-term-mapping). Terms that do not map are searched in all search fields except for Place of Publication, Create Date, Completion Date, Entry Date, MeSH Date, and Modification Date. Terms enclosed in double quotes or including wildcards (*) will be searched in all fields and not processed using automatic term mapping. PubMed ignores [stopwords](#help-stopwords).
 #### Article Identifier [aid]
 Includes article identifiers submitted by journal publishers such as DOI (digital object identifier).
 #### Author [au]
 The format to search for this field is: last name followed by a space and up to the first two initials followed by a space and a suffix abbreviation, if applicable, all without periods or a comma after the last name (e.g., fauci as or o'brien jc jr). Initials and suffixes may be omitted when searching.
 PubMed automatically truncates a search for an author's name to account for varying initials, e.g., o'brien j [au] will retrieve o'brien ja, o'brien jb, o'brien jc jr, as well as o'brien j. To turn off automatic truncation, enclose the author's name in double quotes and tag with [au] in brackets, e.g., "o'brien j" [au] to retrieve just o'brien j.
 Searching by full author name for articles published from 2002 forward is also possible, if available. See [NLM policy on author names](#author-indexing-policy).
 #### Author Identifier [auid]
 The author identifier includes a unique identifier associated with an author, corporate or investigator name, if supplied by a publisher. The field includes the organization authority that established the unique identifier, such as, ORCID, ISNI, VIAF, e.g., orcid 0000-0001-5027-4446 [auid].
 #### Book [book]
 The book search field includes book citations, e.g., genereviews [book].
 Use the following untagged searches to retrieve all book or book chapters, e.g., ataxia AND pmcbookchapter
 - books and chapters: pmcbook
- books: pmcbooktitle
- book chapters: pmcbookchapter
 The above searches capture book records provided by the NCBI Bookshelf database; they exclude a small number of documents from other providers that appear in both PubMed and Bookshelf. For the most comprehensive search of records appearing in both PubMed and Bookshelf, search "pubmed books"[sb].
 #### Comment Correction Type
 The data in these fields are citations to other associated journal publications, e.g., comments or errata. Often these link to the respective citation. Comments/Corrections data can be retrieved by the search term that follows each type:
 - Comment in: hascommentin
- Comment on: hascommenton
- Corrected and republished in: hascorrectedrepublishedin
- Corrected and republished from: hascorrectedrepublishedfrom
- Dataset use reported in: hasassociatedpublication
- Dataset described in: hasassociateddataset
- Erratum in: haserratumin
- Erratum for: haserratumfor
- Expression of concern in: hasexpressionofconcernin
- Expression of concern for: hasexpressionofconcernfor
- Original Report in: hasoriginalreportin
- Republished in: hasrepublishedin
- Republished from: hasrepublishedfrom
- Retracted and republished in: hasretractedandrepublishedin
- Retracted and republished from: hasretractedandrepublishedfrom
- Retraction in: hasretractionin
- Retraction of: hasretractionof
- Summary for patients in: hassummaryforpatientsin
- Update in: hasupdatein
- Update of: hasupdateof
 #### Completion Date [dcom]
 Used by NLM for internal processing. Completon Date is not included in All Fields retrieval; the [dcom] search tag is required.
 #### Conflict of Interest Statement [cois]
 The conflict of interest statement from the published article. Conflict of interest statements are available when supplied by the publisher in the citation data sent to PubMed, or when included in full text articles in PubMed Central (PMC).
 To retrieve all citations that contain conflict of interest statements, use the query hascois.
 #### Corporate Author [cn]
 Corporate author identifies the corporate or collective authorship of an article. Corporate names display exactly as they appear in the journal.
 Note: Citations indexed pre-2000 and some citations indexed in 2000-2001 retain corporate authors at the end of the title field. For comprehensive searches, consider including terms and/or words searched in the title field [ti].
 #### Create Date [crdt]
 The date the citation record was first created in PubMed. Create Date can be helpful when checking PubMed for citations added since the last time a query was run. Create Date is not included in All Fields retrieval; the [crdt] search tag is required.
 #### EC/RN Number [rn]
 EC/RN numbers are assigned by:
 - The Food and Drug Administration (FDA) Substance Registration System for Unique Ingredient Identifiers (UNIIs), e.g., Y92OUS2H9B
- The Enzyme Commission (EC) to designate a particular enzyme, e.g., EC 1.1.1.57
- The Chemical Abstracts Service (CAS) for Registry Numbers, e.g., 2751-14-6
 The EC/RN number search field includes both the Registry Number and the Related Registry Number (available in the NLM MeSH Browser).
 #### Editor [ed]
 The editor search field includes the editors for book or chapter citations.
 #### Entry Date [edat]
 Entry date (EDAT) is used for PubMed processing, such as “Most Recent” sort order (i.e., last in, first out).
 EDAT is typically set within 24 hours of the citation’s availability in PubMed. Exceptions: As of December 15, 2008, citations added to PubMed more than twelve months after the date of publication have the EDAT set to the date of publication, except for book citations. Prior to this, the Entry Date was set to the Publication Date on citations published before September 1997. Entry Date is not included in All Fields retrieval; the [edat] search tag is required.
 Note: Entry Date was called Entrez Date in the legacy PubMed system (retired in 2020).
 #### Filter [filter] [sb]
 Technical tags used by LinkOut, filters include:
 - loall[sb] - citations with LinkOut links in PubMed
- free full text[sb] - citations that include a link to a free full text article
- full text[sb] - citations that include a link to a full text article
 #### First Author Name [1au]
 The first personal author name in a citation.
 #### Full Author Name [fau]
 The full author name for articles published from 2002 forward, if available. Full author searches can be entered in natural or inverted order, e.g., julia s wong or wong julia s.
 #### Full Investigator Name [fir]
 The index for the article's full investigator or collaborator name, if available. Full investigator searches can be entered in natural or inverted order, e.g., harry janes or janes harry.
 #### Grants and Funding [gr]
 The Grants and funding [gr] search field (previously Grant Number) includes grant numbers, contract numbers, or other intramural research identifiers associated with a publication.
 The most common type of funding information associated with a publication in PubMed are grant numbers. Data in the Grants and funding search field can consist of up to four parts:
 1. **Number** contains the grant, contract, intramural project number (or both) that designates financial support by any agency of the United States Public Health Service (US PHS), any institute of the National Institutes of Health, or other organization.
2. **Funder code** contains the 2-letter grant code or funding organization acronym, for example: CA for National Cancer Institute or DDCF for Doris Duke Charitable Foundation. See [Grant Number Information Found in the GR Field in MEDLINE/PubMed *(Archived)*](https://www.nlm.nih.gov/bsd/grant_acronym.html) for the 2-character abbreviations, PHS agency acronyms, and other US and non-US funding organizations.
3. **Agency** includes the acronym or mnemonic in the case of US PHS agencies, or full organization name. As of 2009 this includes the agency's hierarchical structure from lower to higher entity, when known. For example, NCI NIH HHS for National Cancer Institute, National Institutes of Health, U.S. Department of Health and Human Services.
4. **Country** contains the home country of the funding agency, for example: United States.
 Each individual part can be searched using [gr], for example: CA101211[gr], CA[gr], NCI[gr], NIH[gr], or United States[gr].
 This field can also be searched to find articles with intramural support; e.g., "intramural nih"[gr] finds all journal citations authored by intramural NIH staff.
 Completeness of funding information in PubMed will vary by source.
 See [Grants and funding](#grants-and-funding) for more information about data in this field.
 #### Investigator [ir]
 Names of principal investigator(s) or collaborators who contributed to the research. Search names following the author field format, for example: soller b[ir].
 #### ISBN [isbn]
 The ISBN for book or book chapters.
 #### Issue [ip]
 The number of the journal issue in which the article was published.
 #### Journal [ta]
 The journal search field includes the journal title abbreviation, full journal title, or ISSN/eISSN number (e.g., J Biol Chem, Journal of Biological Chemistry, 0021-9258). If a journal title contains special characters, e.g., parentheses, brackets, enter the name without these characters, e.g., enter J Hand Surg [Am] as J Hand Surg Am.
 #### Language [la]
 The language search field includes the language in which the article was published. Note that many non-English articles have English language abstracts. You may search using either the language or the first three characters of most languages, e.g., chi [la] retrieves the same results as chinese [la]. The most notable exception is jpn [la] for Japanese.
 #### Last Author Name [lastau]
 The last personal author name in a citation.
 #### Location ID [lid]
 Location ID includes the DOI or publisher ID that serves the role of pagination to locate an online article.
 #### MeSH Date [mhda]
 The date MeSH terms were added to the citation. MeSH Date is initially set to the [Entry Date [EDAT]](#edat) when the citation is added to PubMed; citations added to PubMed more than twelve months after the date of publication have EDAT and MHDA set to the date of publication.
 MeSH Date is not included in All Fields retrieval; the [mhda] search tag is required.
 Dates must be entered using the format YYYY/MM/DD [mhda], e.g., 2000/03/15 [mhda]. The month and day are optional (e.g., 2000 [mhda] or 2000/03 [mhda]).
 To enter a date range, insert a colon (:) between each date, e.g., 1999:2000 [mhda] or 2000/03:2000/04 [mhda].
 #### MeSH Major Topic [majr]
 A MeSH term that is one of the main topics discussed in the article denoted by an asterisk on the MeSH term or MeSH/Subheading combination, e.g., Cytokines/physiology* See [MeSH Terms [mh]](#mh) below.
 #### MeSH Subheadings [sh]
 [MeSH Subheadings](https://www.nlm.nih.gov/mesh/qualifiers_scopenotes.html) are used with MeSH terms to help describe more completely a particular aspect of a subject. For example, the drug therapy of asthma is displayed as asthma/drug therapy; see MeSH/Subheading Combinations in [MeSH Terms [mh]](#mh) below.
 The MeSH Subheading field allows users to "free float" Subheadings, e.g., hypertension [mh] AND toxicity [sh].
 MeSH Subheadings automatically include the [more specific Subheading terms](https://www.nlm.nih.gov/mesh/subhierarchy.html) under the term in a search. To turn off this automatic feature, use the search syntax [sh:noexp], e.g., therapy [sh:noexp].
 In addition, you can enter the two-letter [MeSH Subheading abbreviations](#mesh-subheadings) rather than spelling out the Subheading, e.g., dh [sh] = diet therapy [sh].
 #### MeSH Terms [mh]
 The [NLM Medical Subject Headings](https://www.nlm.nih.gov/mesh/meshhome.html) controlled vocabulary of biomedical terms that is used to describe the subject of each journal article in MEDLINE. MeSH is updated annually to reflect changes in medicine and medical terminology. MeSH terms are arranged hierarchically by subject categories with more specific terms arranged beneath broader terms. PubMed allows you to view this hierarchy and select terms for searching in the MeSH Database.
 MEDLINE articles are automatically indexed with MeSH terms using a well-refined algorithm. Applying the MeSH vocabulary ensures that articles are uniformly indexed by subject, whatever the author's words. For more information, see [Frequently Asked Questions about Indexing for MEDLINE](https://www.nlm.nih.gov/bsd/indexfaq.html).
 More information about MeSH Terms and Major MeSH Topic search fields:
 - To search the term only as a MeSH term, it must be tagged using the search field, e.g., [mh] for MeSH Terms or [majr] for MeSH Major Topic. A tagged term is checked against the [subject translation table](#automatic-term-mapping), and then mapped to the appropriate MeSH term(s); entry terms tagged with [mh] also map to the appropriate MeSH term(s), increasing the access points to the citations since a searcher may not know the exact preferred term for a MeSH Descriptor. To search for the exact term only and turn off mapping to multiple MeSH terms, enter the tagged MeSH term in double quotes.
- MeSH terms are arranged hierarchically by subject categories with more specific terms arranged beneath broader terms. MeSH terms in PubMed automatically include the more specific MeSH terms in a search. To turn off this automatic feature, use the search syntax [mh:noexp], e.g., neoplasms [mh:noexp].For more detailed information about MeSH vocabulary including the hierarchical structure, please see the [MeSH homepage](https://www.nlm.nih.gov/mesh/meshhome.html).
- MeSH/Subheading Combinations: To directly attach MeSH Subheadings, use the format MeSH Term/Subheading, e.g., neoplasms/diet therapy. You may also use the [two-letter MeSH Subheading abbreviations](#mesh-subheadings), e.g., neoplasms/dh. The [mh] tag is not required, however [majr] may be used, e.g., plants/genetics[majr]. Only one Subheading may be directly attached to a MeSH term. For a MeSH/Subheading combination, PubMed always includes the more specific terms arranged beneath broader terms for the MeSH term and also includes the more specific terms arranged beneath broader [Subheadings](https://www.nlm.nih.gov/mesh/subhierarchy.html). The broader Subheading, or one of its indentions, will be directly attached to the MeSH term or one of its indentions. For example, hypertension/therapy also retrieves hypertension/diet therapy; hypertension/drug therapy; hypertension, malignant/therapy; hypertension, malignant/drug therapy, and so on, as well as hypertension/therapy.
- To turn off the automatic inclusion of the more specific terms, use the syntax [field:noexp], e.g., hypertension [mh:noexp], or hypertension [majr:noexp], or hypertension/therapy [mh:noexp]. The latter example turns off the more specific terms in both parts, searching for only the one Subheading therapy attached directly to only the one MeSH term hypertension.
- If parentheses are embedded in a MeSH term, replace the parentheses with a space and tag with [mh] e.g., enter the MeSH term Benzo(a)pyrene as benzo a pyrene [mh].
- MeSH terms can be selected for searching in the MeSH database and from the advanced search builder index.
 #### Modification Date [lr]
 Modification date is a completed citation’s most recent revision date. Modification Date is not included in All Fields retrieval; the [lr] search tag is required.
 #### NLM Unique ID [jid]
 The NLM ID is the alpha-numeric identifier for the cited journal that was assigned by the NLM Integrated Library System LocatorPlus, e.g., 0375267 [jid].
 #### Other Term [ot]
 The author keyword field (OT field) is searchable with the title/abstract [tiab], text word [tw] and other term [ot] search tags. To retrieve all citations that have keywords, use the query haskeyword. Other term data may display an asterisk to indicate a major concept; however, you cannot search other terms with a major concept tag.
 #### Owner
 The owner search field includes the acronym that identifies the organization that supplied the citation data. Search using owner + the owner acronym, e.g., ownernasa.
 #### Pagination [pg]
 Enter only the first page number that the article appears on. The citation will display the full pagination of the article but this field is searchable using only the first page number.
 #### Personal Name as Subject [ps]
 Use this search field tag to limit retrieval to where the name is the subject of the article, e.g., varmus h[ps]. Search for personal names as subject using the author field format, e.g., varmus h[ps].
 #### Pharmacological Action [pa]
 Substances known to have a particular pharmacologic action. Each pharmacologic action term index is created with the drug/substance terms known to have that effect. This includes both MeSH terms and terms for Supplementary Concept Records.
 #### Place of Publication [pl]
 Indicates the cited journal's country of publication. Geographic place of publication regions are not searchable. In order to retrieve records for all countries in a region (e.g., North America) it is necessary to OR together the countries of interest. Note: This field is not included in all fields or text word retrieval.
 #### PMCID and MID
 Search for PMC or NIH manuscript identifiers using the appropriate prefix followed by the ID number, e.g., PMC2600426. To retrieve all NIH manuscript citations, use the query hasnihmsid.
 #### PMID [pmid]
 To search for a PubMed Identifier (PMID), enter the ID with or without the search field tag [pmid]. You can search for several PMIDs by entering each number in the search box separated by a space (e.g., 17170002 16381840); PubMed will OR the PMIDs together. To search PMIDs in combination with other terms, the search field tag is required: lipman AND 16381840[pmid].
 PMIDs do not change over time or during processing and are never reused.
 #### Publication Date [dp]
 Publication date is the date that the article was published. The search field tags [dp] and [pdat] may be used interchangeably for publication date searching.
 Dates or date ranges must be searched using the format yyyy/mm/dd [dp], e.g., 1998/03/06 [dp]. The month and day are optional (e.g., 1998 [dp] or 1998/03 [dp]).
 To enter a date range search, insert a colon (:) between each date, e.g., 1996:1998 [dp] or 1998/01:1998/04 [dp].
 Use the following format to search X days, months or years immediately preceding today’s date where X = numeric value:
 - "last X days"[dp]
- "last X months"[dp]
- "last X year"[dp]
 More information about publication dates:
 - The time between an article’s publication and the citation’s availability in PubMed varies depending on when the publisher deposits the citation to PubMed. Because of this, searching with Create Date [crdt] (the date a citation was created in PubMed) is often more comprehensive than Publication Date [dp] when checking PubMed on a regular basis for new citations.
- Journals vary in the way the publication date appears on an issue. Some journals include just the year, whereas others include the year plus month or year plus month plus day. And, some journals use the year and season (e.g., Winter 1997). The publication date in the citation is recorded as it appears in the journal.
- Publication dates without a month are set to January, multiple months (e.g., Oct-Dec) are set to the first month, and dates without a day are set to the first day of the month. Dates with a season are set as: winter = January, spring = April, summer = July and fall = October.
- If an article is published electronically and in print on different dates both dates are searchable and may be included on the citation prefaced with an Epub or Print label. Note: The electronic date is not included in a Publication Date [DP] search if the electronic date is later than the print date.
- To search for electronic dates only use the search tag [EPDAT], for print dates only tag with [PPDAT].
- Most journals now publish articles online on a continuous basis, as soon as they are ready for publication (after peer review and editing, etc.) instead of, or in addition to, publishing collections of articles as an "issue" on a periodic basis. When a journal deposits a citation for an "online first" article in PubMed, NLM appends the note "[Online ahead of print]" to the online publication date. The citation is updated, and the ahead of print notation removed, when the article is included in a journal issue. The lag between the "online first" and "issue" publication dates may be days, weeks, months, or more than a year. In many cases, depending on the journal, the online first version is considered to be the version of record. The "[Online ahead of print]" note in PubMed should not be taken to mean that the cited article is not the version of record.
- Bookshelf citation publication dates are generated from the book’s original publication date. Publication date searches include Contribution Date (CTDT) and Date Revised (DRDT) in addition to Publication Date (DP) for Book citations.
 #### Publication Type [pt]
 Describes the material presented in the article (e.g., Review, Clinical Trial, Retracted Publication, Letter). Citations may include multiple Publication Types. Use the search tag [pt] with any PubMed [Publication Type](#publication-types), e.g., review[pt].
 Publication Types are arranged hierarchically with more specific terms arranged beneath broader terms, and publication types automatically include the more specific publication types in a search.
 - To turn off this automatic feature, use the search syntax [pt:noexp], e.g., review [pt:noexp].
- To view the publication type hierarchy, see the [Publication Formats: MeSH Tree Structures](https://meshb.nlm.nih.gov/record/ui?ui=D052180) in the MeSH Browser.
 #### Publisher [pubn]
 Includes publisher names for Bookshelf citations.
 #### Secondary Source ID [si]
 The SI field identifies secondary source databanks and accession numbers, e.g., GenBank, GEO, PubChem, [ClinicalTrials.gov](https://ClinicalTrials.gov), ISRCTN. The field is composed of the source followed by a slash followed by an accession number and can be searched with one or both components, e.g., genbank [si], AF001892 [si], genbank/AF001892 [si]. To retrieve all citations with an SI value, search hasdatabanklist.
 #### Subset [sb]
 The subset field is a method of restricting retrieval by subject, citation status and journal category, with the search tag [SB]. See also [filters](#help-filters) and [Find related resources using LinkOut.](#using-linkout)
 #### Supplementary Concept [nm]
 Includes chemical, protocol, disease or organism terms. Synonyms to the supplementary concepts will automatically map when tagged with [nm]. This field was implemented in mid-1980; however, many chemical names are searchable as MeSH terms before that date.
 #### Text Words [tw]
 Includes all words and numbers in the title, abstract, other abstract, MeSH terms, MeSH Subheadings, Publication Types, Substance Names, Personal Name as Subject, Corporate Author, Secondary Source, Comment/Correction Notes, and Other Terms (see [Other Term [OT]](#ot) above) typically non-MeSH subject terms (keywords), including NASA Space Flight Mission, assigned by an organization other than NLM.
 #### Title [ti]
 Words and numbers included in the title of a citation, as well as the collection title for book citations.
 #### Title/Abstract [tiab]
 Words and numbers included in a citation's title, collection title, abstract, other abstract and author keywords ( [Other Term [ot]](#ot) field). English language abstracts are taken directly from the published article. If an article does not have a published abstract, NLM does not create one.
 #### Transliterated Title [tt]
 Words and numbers in title originally published in a non-English language, in that language. Non-Roman alphabet language titles may be transliterated. Transliterated title is not included in [Text Word [TW]](#tw) retrieval.
 #### Volume [vi]
 The number of the journal volume in which an article is published.
 ### NLM author indexing policy
 NLM author indexing policy is as follows:
 - 1966 - 1984: MEDLINE did not limit the number of authors.
- 1984 - 1995: The NLM limited the number of authors to 10, with "et al" as the eleventh occurrence.
- 1996 - 1999: The NLM increased the limit from 10 to 25. If there were more than 25 authors, the first 24 were listed, the last author was used as the 25th, and the twenty-sixth and beyond became "et al."
- 2000 - Present: MEDLINE does not limit the number of authors.
 More information:
 - Beginning in mid-2005, the policy restrictions on number of author names in past years were lifted so that on an individual basis, a citation may be edited to include all author names in the published article, regardless of the limitation in effect when the citation was created.
- Effective with 1992 date of publication, letters are indexed individually with authors rather than as an anonymous group.
- Until 1990, NLM transliterated up to five authors' Cyrillic or Japanese names to the Roman alphabet.
- Between 1990 and 2016, the first ten Cyrillic or Japanese names are transliterated. Chinese ideograms were not transliterated by NLM, but if transliterations of the authors names are available in the journal article or table of contents, they were included in the citation, even if that includes only one author in a multi-author article.
- Beginning in 2016, author names are published in Roman characters in all MEDLINE journals, and NLM no longer transliterates Cyrillic or Japanese names. All author names are included as published.
 ### Error messages
 #### System error messages
 Please provide your IT staff with the technical [browser advice for NCBI web pages](https://www.ncbi.nlm.nih.gov/guide/browsers/) to ensure your browser, firewall, and servers are enabled for JavaScript, cookies, pop-ups, and HTTP 1.1. Antivirus software may affect page caching which can result in unexpected page expired messages. Also, nlm.nih.gov should be added as a browser exception and be considered a trusted site by your system and network. You may have to delete your browser's cache (temporary files) before trying to access PubMed again.
 #### Typographical errors
 Please contact the journal publisher directly to report an error and initiate a correction to PubMed citations for content other than MeSH.
 To report a MeSH error in a PubMed citation, please contact the [NLM Help Desk](https://support.nlm.nih.gov/support/create-case/?category=pubmed) and include the PMID number (e.g., PMID: 12345678), and an indication of the incorrect and correct information.
 NLM provides data to vendors around the world. Other products and services will not necessarily immediately reflect corrections made to PubMed records. If you search through a vendor's system, please contact your vendor about their maintenance schedules.
 ### Cookies
 A "cookie" is information stored by a web site server on your computer. See the [NLM Privacy Policy](https://www.nlm.nih.gov/privacy.html) for additional information.
 In the case of PubMed, cookies store information about your interactions that may be needed later to perform a function. To use these interactive features you need to enable cookies on your computer. Consult your browser's help for information on enabling cookies.
 If you have problems using cookie-dependent features of PubMed, even after enabling cookies, possible reasons may include:
 - Cookies are blocked by your provider or institution. Check with your Internet provider and/or the system administrator at your institution to see if cookies can be accepted. Even if you have them enabled in your web browser, if they are blocked by your provider or institution (e.g., by a firewall, proxy server, etc.), cookie-dependent features of PubMed won't work.
- Your computer's date and time settings are incorrect. Check your computer's time settings to ensure that they are correct.
 ### MeSH Subheadings
 See the MeSH Subheadings table below and [scope notes and allowable categories](https://www.nlm.nih.gov/mesh/qualifiers_scopenotes.html) on the NLM website.
 Abbreviation MeSH Subheading Abbreviation MeSH Subheading AB Abnormalities IR Innervation AD Administration and Dosage IS Instrumentation AE Adverse Effects IP Isolation and Purification AG Agonists LJ Legislation and Jurisprudence AA Analogs and Derivatives ME Metabolism AN Analysis MT Methods AH Anatomy and Histology MI Microbiology AI Antagonists and Inhibitors MO Mortality BI Biosynthesis NU Nursing BS Blood Supply OG Organization and Administration BL Blood PS Parasitology CF Cerebrospinal Fluid PY Pathogenicity CS Chemical Synthesis PA Pathology CI Chemically Induced PK Pharmacokinetics CH Chemistry PD Pharmacology CL Classification PH Physiology CO Complications PP Physiopathology CN Congenital PO Poisoning CY Cytology PC Prevention and Control DF Deficiency PX Psychology DI Diagnosis RE Radiation Effects DH Diet Therapy RT Radiotherapy DG Diagnostic Imaging RH Rehabilitation DE Drug Effects SC Secondary DT Drug Therapy ST Standards EC Economics SN Statistics and Numerical Data ED Education SD Supply and Distribution EM Embryology SU Surgery EN Enzymology TU Therapeutic Use EP Epidemiology TH Therapy ES Ethics TO Toxicity EH Ethnology TM Transmission ET Etiology TR Transplantation GE Genetics TD Trends GD Growth and Development UL Ultrastructure HI History UR Urine IM Immunology VE Veterinary IN Injuries VI Virology ### Stopwords
 Stopwords A a, about, again, all, almost, also, although, always, among, an, and, another, any, are, as, at B be, because, been, before, being, between, both, but, by C can, could D did, do, does, done, due, during E each, either, enough, especially, etc F for, found, from, further H had, has, have, having, here, how, however I i, if, in, into, is, it, its, itself J just K kg, km M made, mainly, make, may, mg, might, ml, mm, most, mostly, must N nearly, neither, no, nor O obtained, of, often, on, our, overall P perhaps, pmid Q quite R rather, really, regarding S seem, seen, several, should, show, showed, shown, shows, significantly, since, so, some, such T than, that, the, their, theirs, them, then, there, therefore, these, they, this, those, through, thus, to U upon V various, very W was, we, were, what, when, which, while, with, within, without, would ### PubMed character conversions
 Certain characters have special meaning in searches, others are converted to spaces.
 Searches that include the following characters are translated as follows:
 - parentheses ( ) - used to create Boolean nesting
- square brackets [ ] - search field tag qualification
- ampersand & - Boolean operator AND
- pipe | - Boolean operator OR
- forward slash / - MeSH/Subheading combinations
- colon : - designates a range operation
- double quotes " - used to force a phrase search
- pound sign # - designates a History search statement when immediately followed by a number, e.g., #1 AND cat
- asterisk * - wildcard symbol, e.g., toxicol*
 Characters converted to spaces in search queries:
 - exclamation mark !
- pound sign #
- dollar sign $
- percentage sign %
- asterisk * (if it cannot be used in a wildcard search, for example, when a term is too short)
- plus symbol +
- minus symbol -
- period .
- comma ,
- semi-colon ;
- angle brackets < >
- equal sign =
- question mark ?
- backslash \
- caret ^
- underscore _
- curly brackets { }
- approximately ~
- single quotes '
 Some characters have special meaning in [MeSH fields:](https://www.nlm.nih.gov/mesh/meshhome.html)
 - forward slash /
- hypens -
- comma ,
 ### Publication Types
 Publication types found in PubMed are listed below. For more information, see [Publication Type [PT]](#pt) and [MeSH Publication Types with Scope Notes](https://www.nlm.nih.gov/mesh/pubtypes.html); however, please note that not all MeSH publication types are used in PubMed, and some that were previously used have been discontinued.
 Discontinued publication types are marked in the list below with an asterisk (*). These publication types remain on previously indexed citations and are searchable in PubMed using the [[pt] search field tag](#pt); however, they are not applied to new citations and do not appear in the list of [article type filters](#filters-article-type) available from the filter sidebar.
 - Adaptive Clinical Trial
- Address
- Autobiography *
- Bibliography *
- Biography
- Case Reports
- Clinical Conference *
- Clinical Study
- Clinical Trial
- Clinical Trial, Phase I
- Clinical Trial, Phase II
- Clinical Trial, Phase III
- Clinical Trial, Phase IV
- Clinical Trial Protocol
- Clinical Trial, Veterinary
- Collected Work
- Comment
- Comparative Study
- Conference Proceedings
- Consensus Statement
- Consensus Development Conference, NIH *
- Controlled Clinical Trial
- Corrected and Republished Article
- Dataset
- Dictionary *
- Directory *
- Duplicate Publication
- Editorial
- Electronic Supplementary Materials
- English Abstract
- Equivalence Trial
- Evaluation Study
- Evidence Synthesis
- Expression of Concern
- Festschrift
- Government Publication *
- Guideline
- Historical Article
- Interactive Tutorial *
- Interview
- Introductory Journal Article
- Journal Article (Default value when no more descriptive PT is provided or assigned)
- Lecture
- Legal Case *
- Legislation *
- Letter
- Meta-Analysis
- Multicenter Study
- Network Meta-Analysis
- News
- Newspaper Article *
- Observational Study
- Observational Study, Veterinary
- Overall *
- Patient Education Handout
- Periodical Index *
- Personal Narrative
- Portrait *
- Practice Guideline
- Pragmatic Clinical Trial
- Preprint
- Published Erratum
- Randomized Controlled Trial
- Randomized Controlled Trial, Veterinary
- Research Support, American Recovery and Reinvestment Act
- Research Support, N.I.H., Extramural
- Research Support, N.I.H., Intramural
- Research Support, Non-U.S. Gov't
- Research Support, U.S. Gov't, Non-P.H.S.
- Research Support, U.S. Gov't, P.H.S.
- Retracted Publication
- Retraction Notice
- Review
- Scientific Integrity Review *
- Scoping Review
- Seminal Article *
- Systematic Review
- Technical Report *
- Twin Study
- Validation Study
- Video-Audio Media
- Webcast
 ### Status Subsets
 How to Search Citation Status publisher[sb] NOT pubstatusnihms NOT pubstatuspmcsd NOT pmcbook Citations recently added to PubMed via electronic submission from a publisher, and are soon to proceed to the next stage, PubMed - in process (see below). Also for citations received before late 2003 if they are from journals not indexed for MEDLINE, or from a journal that was accepted for MEDLINE after the citations' publication date. These citations bibliographic data have not been reviewed. inprocess[sb] MeSH terms will be assigned if the subject of the article is within the scope of MEDLINE. medline[sb] Citations that have been indexed with MeSH terms, Publication Types, Substance Names, etc. pubstatusnihms AND publisher[sb] Author manuscripts submitted to PMC that fall under the NIH Public Access Policy. pubstatuspmcsd AND publisher[sb] Records for selective deposit articles in PMC. These are articles published in non-MEDLINE journals where the publisher has chosen to deposit in PMC only those articles that fall under the NIH Public Access Policy. pmcbook Book and book chapter citations available on the NCBI Bookshelf. pubmednotmedline[sb] Citations that will not receive MEDLINE indexing because they are for articles in non-MEDLINE journals, or they are for articles in MEDLINE journals but the articles are out of scope, or they are from issues published prior to the date the journal was selected for indexing, or citations to articles from journals that deposit their full text articles in PMC but have not yet been recommended for indexing in MEDLINE. ### Filter search strategies
 - [Publication date](#filter-strategy-pubdate)
- [Text availability](#filter-strategy-text-availability)
- [Article attribute](#filter-strategy-article-attribute)
- [Article type](#filter-strategy-article-type)
- [Article language](#filter-strategy-language)
- [Species](#filter-strategy-species)
- [Sex](#filter-strategy-sex)
- [Age](#filter-strategy-age)
- [Other](#filter-strategy-other)
 #### Publication date
 Filter name PubMed equivalent 1 year "last year"[dp] 5 years "last 5 years"[dp] 10 years "last 10 years"[dp] Custom range See: [Searching for a date range](#date-search-range) #### Text availability
 Filter name PubMed equivalent Abstract hasabstract Free full text free full text[sb] Full text full text[sb] #### Article attribute
 Filter name PubMed equivalent Associated data data[sb] #### Article type
 Most article type filters use the article type name with the [publication type [pt]](#pt) search field tag; for example, "multicenter study"[pt]. Exceptions are noted in the table below.
 Filter name PubMed equivalent Books and Documents "pubmed books"[sb] Systematic Review Uses a [search strategy](https://www.nlm.nih.gov/bsd/pubmed_subsets/sysreviews_strategy.html) that includes the publication type [pt] Retraction Notice hasretractionof OR retraction notice[pt] Retracted Publication hasretractionin OR retracted publication[pt] #### Article language
 The [article language filters](#filters-language) use the language name with the [language [la]](#la) search field tag; for example, esperanto[la].
 #### Species
 Filter name PubMed equivalent Humans humans[mh] Other animals "animals"[mh:noexp] #### Sex
 Filter name PubMed equivalent Female female[mh] Male male[mh] #### Age
 Filter name PubMed equivalent Child: birth-18 years (infant[mh] OR child[mh] OR adolescent[mh]) Newborn: birth-1 month infant, newborn[mh] Infant: birth-23 months infant[mh] Infant: 1-23 months infant[mh:noexp] Preschool Child: 2-5 years child, preschool[mh] Child: 6-12 years child[mh:noexp] Adolescent: 13-18 years adolescent[mh] Adult: 19+ years adult[mh] Young Adult: 19-24 years "young adult"[mh] Adult: 19-44 years adult[mh:noexp] Middle Aged + Aged: 45+ years (middle aged[mh] OR aged[mh]) Middle Aged: 45-64 years middle aged[mh] Aged: 65+ years aged[mh] 80 and over: 80+ years aged, 80 and over[mh] #### Other
 See [Other filters and more subsets](#filters-other).
 ### Clinical Queries filters
 - [COVID-19 Articles](#covid19-article-filters)
- [Clinical Study Categories](#clinical-study-category-filters)
- [Medical Genetics](#medical-genetics-filters)
 #### COVID-19 article filters
 The COVID-19 article filters limit retrieval to citations about the 2019 novel coronavirus; these filters may evolve over time.
 Category Filter name PubMed equivalent General LitCGeneral ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) Mechanism LitCMechanism ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("metabolic"[All Fields] OR "metabolical"[All Fields] OR "metabolically"[All Fields] OR "metabolics"[All Fields] OR "metabolism"[MeSH Terms] OR "metabolism"[All Fields] OR "metabolisms"[All Fields] OR "metabolism"[MeSH Subheading] OR "metabolic networks and pathways"[MeSH Terms] OR ("metabolic"[All Fields] AND "networks"[All Fields] AND "pathways"[All Fields]) OR "metabolic networks and pathways"[All Fields] OR "metabolities"[All Fields] OR "metabolization"[All Fields] OR "metabolize"[All Fields] OR "metabolized"[All Fields] OR "metabolizer"[All Fields] OR "metabolizers"[All Fields] OR "metabolizes"[All Fields] OR "metabolizing"[All Fields] OR "virology"[MeSH Subheading] OR ("mechanism"[All Fields] OR "mechanisms"[All Fields]) OR ("etiology"[MeSH Subheading] OR "etiology"[All Fields] OR "pathogenesis"[All Fields]) OR "pathologic process*"[All Fields]) Transmission LitCTransmission ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("transmission"[Text Word] OR "transmission"[MeSH Subheading] OR "replication"[Text Word] OR "disease transmission, infectious"[MeSH Terms]) Diagnosis LitCDiagnosis ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("diagnos*"[All Fields] OR "detect*"[All Fields] OR "diagnosis"[MeSH Terms] OR "diagnostic equipment"[MeSH Terms] OR "diagnostic errors"[MeSH Terms] OR "diagnostic imaging"[MeSH Terms] OR "diagnostic services"[MeSH Terms] OR "diagnosis, differential"[MeSH Terms] OR "diagnosis"[MeSH Subheading]) Treatment LitCTreatment ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("therapeutics"[MeSH Terms] OR "therapeutics"[All Fields] OR "treatments"[All Fields] OR "therapy"[MeSH Subheading] OR "therapy"[All Fields] OR "treatment"[All Fields] OR "treatment s"[All Fields] OR "treat*"[All Fields] OR ("clinical trial"[Publication Type] OR "clinical trials as topic"[MeSH Terms] OR "clinical trials"[All Fields]) OR ("clinical trial"[Publication Type] OR "clinical trials as topic"[MeSH Terms] OR "clinical trial"[All Fields]) OR ("randomized controlled trial"[Publication Type] OR "randomized controlled trials as topic"[MeSH Terms] OR "randomized controlled trial"[All Fields] OR "randomised controlled trial"[All Fields]) OR ("randomized controlled trial"[Publication Type] OR "randomized controlled trials as topic"[MeSH Terms] OR "randomized controlled trials"[All Fields] OR "randomised controlled trials"[All Fields]) OR ("therapeutics"[MeSH Terms] OR "therapeutics"[All Fields] OR "therapies"[All Fields] OR "therapy"[MeSH Subheading] OR "therapy"[All Fields] OR "therapy s"[All Fields] OR "therapys"[All Fields]) OR ("therapeutical"[All Fields] OR "therapeutically"[All Fields] OR "therapeuticals"[All Fields] OR "therapeutics"[MeSH Terms] OR "therapeutics"[All Fields] OR "therapeutic"[All Fields])) Prevention LitCPrevention ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("transmission*"[All Fields] OR "prevent*"[All Fields] OR "intervent*"[All Fields] OR ("prognosis"[MeSH Terms] OR "prognosis"[All Fields] OR "prognoses"[All Fields]) OR "treatment outcome"[All Fields] OR "prevention and control"[MeSH Subheading] OR ("therapeutical"[All Fields] OR "therapeutically"[All Fields] OR "therapeuticals"[All Fields] OR "therapeutics"[MeSH Terms] OR "therapeutics"[All Fields] OR "therapeutic"[All Fields]) OR "therapeutic*"[All Fields]) Case Report LitCCaseReport ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("case report*"[All Fields] OR ("case reports"[Publication Type] OR "case reports"[All Fields]) OR "report a case"[All Fields] OR ("report*"[All Fields] AND ("ambulatory care facilities"[MeSH Terms] OR ("ambulatory"[All Fields] AND "care"[All Fields] AND "facilities"[All Fields]) OR "ambulatory care facilities"[All Fields] OR "clinic"[All Fields] OR "clinic s"[All Fields] OR "clinical"[All Fields] OR "clinically"[All Fields] OR "clinicals"[All Fields] OR "clinics"[All Fields] OR "patient*"[All Fields])) OR "reported case"[All Fields] OR "clinical presentation*"[All Fields] OR "patient management"[All Fields] OR "infected patient*"[All Fields]) Forecasting LitCForecasting ("COVID-19" OR "COVID19" OR "COVID-19"[MeSH Terms] OR "COVID-19 Vaccines" OR "COVID-19 Vaccines"[MeSH Terms] OR "COVID-19 serotherapy" OR "COVID-19 serotherapy"[ MeSH Terms] OR "COVID-19 Nucleic Acid Testing" OR "covid-19 nucleic acid testing"[MeSH Terms] OR "COVID-19 Serological Testing" OR "covid-19 serological testing"[MeSH Terms] OR "COVID-19 Testing" OR "covid-19 testing"[MeSH Terms] OR "SARS-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "sars-cov-2"[MeSH Terms] OR "Severe Acute Respiratory Syndrome Coronavirus 2" OR "2019 NCOV" OR (("coronavirus"[MeSH Terms] OR "coronavirus" OR "COV" OR "NCOV") AND 2019/11/01[PDAT] : 3000/12/31[PDAT])) AND ("forecast*"[All Fields] OR ("forecasted"[All Fields] OR "forecaster"[All Fields] OR "forecasters"[All Fields] OR "forecasting"[MeSH Terms] OR "forecasting"[All Fields] OR "forecast"[All Fields] OR "forecasts"[All Fields] OR "trends"[MeSH Subheading] OR "trends"[All Fields]) OR "trend*"[All Fields] OR "prediction*"[All Fields]) Long COVID LitCLongCOVID "COVID-19 sequela*" OR (("COVID-19" OR "COVID19" OR "Sars-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "2019 Novel Coronavirus" OR "2019-nCoV" OR "Coronavirus Disease 2019" OR "Coronavirus Disease-19" OR "SARS Coronavirus 2" OR "Severe Acute Respiratory Syndrome Coronavirus 2") AND sequela*) OR "post acute sequelae of Sars-CoV-2" OR ("PASC" AND ("COVID-19" OR "COVID19" OR "Sars-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "2019 Novel Coronavirus" OR "2019-nCoV" OR "Coronavirus Disease 2019" OR "Coronavirus Disease-19" OR "SARS Coronavirus 2" OR "Severe Acute Respiratory Syndrome Coronavirus 2")) OR "post acute sequelae of COVID" OR (("post-intensive care syndrome" OR "postintensive care syndrome") AND ("COVID-19" OR "COVID19" OR "Sars-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "2019 Novel Coronavirus" OR "2019-nCoV" OR "Coronavirus Disease 2019" OR "Coronavirus Disease-19" OR "SARS Coronavirus 2" OR "Severe Acute Respiratory Syndrome Coronavirus 2")) OR "post COVID condition*" OR ("PCC" AND ("COVID-19" OR "COVID19" OR "Sars-CoV-2" OR "SARSCoV2" OR "SARSCoV-2" OR "SARS-CoV2" OR "2019 Novel Coronavirus" OR "2019-nCoV" OR "Coronavirus Disease 2019" OR "Coronavirus Disease-19" OR "SARS Coronavirus 2" OR "Severe Acute Respiratory Syndrome Coronavirus 2")) OR "convalescent COVID-19" OR "long haul COVID" OR "COVID long haul*" OR "long COVID" OR "long term COVID" OR "COVID-19 survivor*" OR "post COVID-19 symptom*" OR "chronic COVID syndrome" OR "post COVID syndrome" OR "post COVID-19 neurological syndrome" OR "post acute COVID-19" OR "post-acute COVID-19 syndrome"[MeSH Terms] OR "COVID-19 post-intensive care syndrome"[Supplementary Concept] #### Clinical Study Categories
 The Clinical Study Categories search filters are based on the work of [Haynes RB et al.](#clinical-study-categories-bibliography)
 Category Optimized for Sensitive/ Specific PubMed equivalent Therapy sensitive/broad 99%/70% ((clinical[Title/Abstract] AND trial[Title/Abstract]) OR clinical trials as topic[MeSH Terms] OR clinical trial[Publication Type] OR random*[Title/Abstract] OR random allocation[MeSH Terms] OR therapeutic use[MeSH Subheading]) specific/narrow 93%/97% (randomized controlled trial[Publication Type] OR (randomized[Title/Abstract] AND controlled[Title/Abstract] AND trial[Title/Abstract])) Diagnosis sensitive/broad 98%/74% (sensitiv*[Title/Abstract] OR sensitivity and specificity[MeSH Terms] OR diagnose[Title/Abstract] OR diagnosed[Title/Abstract] OR diagnoses[Title/Abstract] OR diagnosing[Title/Abstract] OR diagnosis[Title/Abstract] OR diagnostic[Title/Abstract] OR diagnosis[MeSH:noexp] OR (diagnostic equipment[MeSH:noexp] OR diagnostic errors[MeSH:noexp] OR diagnostic imaging[MeSH:noexp] OR diagnostic services[MeSH:noexp]) OR diagnosis, differential[MeSH:noexp] OR diagnosis[Subheading:noexp]) specific/narrow 64%/98% (specificity[Title/Abstract]) Etiology sensitive/broad 93%/63% (risk*[Title/Abstract] OR risk*[MeSH:noexp] OR (risk adjustment[MeSH:noexp] OR risk assessment[MeSH:noexp] OR risk factors[MeSH:noexp] OR risk management[MeSH:noexp] OR risk taking[MeSH:noexp]) OR cohort studies[MeSH Terms] OR group[Text Word] OR groups[Text Word] OR grouped [Text Word]) specific/narrow 51%/95% ((relative[Title/Abstract] AND risk*[Title/Abstract]) OR (relative risk[Text Word]) OR risks[Text Word] OR cohort studies[MeSH:noexp] OR (cohort[Title/Abstract] AND study[Title/Abstract]) OR (cohort[Title/Abstract] AND studies[Title/Abstract])) Prognosis sensitive/broad 90%/80% (incidence[MeSH:noexp] OR mortality[MeSH Terms] OR follow up studies[MeSH:noexp] OR prognos*[Text Word] OR predict*[Text Word] OR course*[Text Word]) specific/narrow 52%/94% (prognos*[Title/Abstract] OR (first[Title/Abstract] AND episode[Title/Abstract]) OR cohort[Title/Abstract]) Clinical Prediction Guides sensitive/broad 96%/79% (predict*[Title/Abstract] OR predictive value of tests[MeSH Terms] OR score[Title/Abstract] OR scores[Title/Abstract] OR scoring system[Title/Abstract] OR scoring systems[Title/Abstract] OR observ*[Title/Abstract] OR observer variation[MeSH Terms]) specific/narrow 54%/99% (validation[Title/Abstract] OR validate[Title/Abstract]) ##### Clinical Study Categories bibliography
 The [Clinical Queries search strategies](#clinical-study-category-filters) have been updated based on new evidence from Haynes et al. The current strategies have better performance than their [predecessors](https://pubmed.ncbi.nlm.nih.gov/7850570/). Details of methods appear in the references below.
 Revised December 2011
 - Wilczynski NL, McKibbon KA, Haynes RB. Sensitive Clinical Queries retrieved relevant systematic reviews as well as primary studies: an analytic survey. J Clin Epidemiol. 2011 Dec;64(12):1341-9. doi: 10.1016/j.jclinepi.2011.04.007. Epub 2011 Jul 19. [PMID: 21775104](https://pubmed.ncbi.nlm.nih.gov/21775104/).
- Lokker C, Haynes RB, Wilczynski NL, McKibbon KA, Walter SD. Retrieval of diagnostic and treatment studies for clinical use through PubMed and PubMed's Clinical Queries filters. J Am Med Inform Assoc. 2011 Sep-Oct;18(5):652-9. doi: 10.1136/amiajnl-2011-000233. Epub 2011 Jun 15. [PMID: 21680559](https://pubmed.ncbi.nlm.nih.gov/21680559/); [PMCID: PMC3168323](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3168323/).
- Wilczynski NL, Haynes RB; QI Hedges Team. Optimal search filters for detecting quality improvement studies in Medline. Qual Saf Health Care. 2010 Dec;19(6):e31. doi: 10.1136/qshc.2010.042432. Epub 2010 Jul 29. [PMID: 20671080](https://pubmed.ncbi.nlm.nih.gov/20671080/).
- Kastner M, Wilczynski NL, McKibbon AK, Garg AX, Haynes RB. Diagnostic test systematic reviews: bibliographic search filters ("Clinical Queries") for diagnostic accuracy studies perform well. J Clin Epidemiol. 2009 Sep;62(9):974-81. doi: 10.1016/j.jclinepi.2008.11.006. Epub 2009 Feb 20. PMID: 19230607; [PMCID: PMC2737707](https://pubmed.ncbi.nlm.nih.gov/19230607/).
- Wilczynski NL, Haynes RB. Response to Corrao et al.: Improving efficacy of PubMed clinical queries for retrieving scientifically strong studies on treatment. J Am Med Inform Assoc. 2007 Mar-Apr;14(2):247-8. Epub 2007 Jan 9. [PMID: 17213490](https://pubmed.ncbi.nlm.nih.gov/17213490/); [PMCID: PMC2213472](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2213472/).
- Wilczynski NL, McKibbon KA, Haynes RB. Response to Glanville et al.: How to identify randomized controlled trials in MEDLINE: ten years on. J Med Libr Assoc. 2007 Apr;95(2):117-8; author reply 119-20. [PMID: 17443240](https://pubmed.ncbi.nlm.nih.gov/17443240/); [PMCID: PMC1852612](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1852612/).
- Wilczynski NL, Morgan D, Haynes RB; Hedges Team. An overview of the design and methods for retrieving high-quality studies for clinical care. BMC Med Inform Decis Mak. 2005 Jun 21;5:20. doi: 10.1186/1472-6947-5-20. [PMID: 15969765](https://pubmed.ncbi.nlm.nih.gov/15969765/); [PMCID: PMC1183213](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1183213/).
- Haynes RB, McKibbon KA, Wilczynski NL, Walter SD, Werre SR; Hedges Team. Optimal search strategies for retrieving scientifically strong studies of treatment from Medline: analytical survey. BMJ. 2005 May 21;330(7501):1179. doi: 10.1136/bmj.38446.498542.8F. Epub 2005 May 13. [PMID: 15894554](https://pubmed.ncbi.nlm.nih.gov/15894554/); [PMCID: PMC558012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC558012/).
- Montori VM, Wilczynski NL, Morgan D, Haynes RB; Hedges Team. Optimal search strategies for retrieving systematic reviews from Medline: analytical survey. BMJ. 2005 Jan 8;330(7482):68. doi: 10.1136/bmj.38336.804167.47. Epub 2004 Dec 24. [PMID: 15619601](https://pubmed.ncbi.nlm.nih.gov/15619601/); [PMCID: PMC543864](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC543864/).
- Wilczynski NL, Haynes RB, Lavis JN, Ramkissoonsingh R, Arnold-Oatley AE; HSR Hedges team. Optimal search strategies for detecting health services research studies in MEDLINE. CMAJ. 2004 Nov 9;171(10):1179-85. doi: 10.1503/cmaj.1040512. [PMID: 15534310](https://pubmed.ncbi.nlm.nih.gov/15534310/); [PMCID: PMC524948](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC524948/).
- Wilczynski NL, Haynes RB; Hedges Team. Developing optimal search strategies for detecting clinically sound prognostic studies in MEDLINE: an analytic survey. BMC Med. 2004 Jun 9;2:23. doi: 10.1186/1741-7015-2-23. [PMID: 15189561](https://pubmed.ncbi.nlm.nih.gov/15189561/); [PMCID: PMC441418](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC441418/).
- Haynes RB, Wilczynski NL. Optimal search strategies for retrieving scientifically strong studies of diagnosis from Medline: analytical survey. BMJ. 2004 May 1;328(7447):1040. doi: 10.1136/bmj.38068.557998.EE. Epub 2004 Apr 8. [PMID: 15073027](https://pubmed.ncbi.nlm.nih.gov/15073027/); [PMCID: PMC403841](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC403841/).
- Bhandari M, Montori VM, Devereaux PJ, Wilczynski NL, Morgan D, Haynes RB; Hedges Team. Doubling the impact: publication of systematic review articles in orthopaedic journals. J Bone Joint Surg Am. 2004 May;86(5):1012-6. [PMID: 15118046](https://pubmed.ncbi.nlm.nih.gov/15118046/).
- Wong SS, Wilczynski NL, Haynes RB; Hedges Team. Developing optimal search strategies for detecting clinically relevant qualitative studies in MEDLINE. Stud Health Technol Inform. 2004;107(Pt 1):311-6. [PMID: 15360825](https://pubmed.ncbi.nlm.nih.gov/15360825/).
- Montori VM, Wilczynski NL, Morgan D, Haynes RB; Hedges Team. Systematic reviews: a cross-sectional study of location and citation counts. BMC Med. 2003 Nov 24;1:2. doi: 10.1186/1741-7015-1-2. [PMID: 14633274](https://pubmed.ncbi.nlm.nih.gov/14633274/); [PMCID: PMC281591](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC281591/).
- Wong SS, Wilczynski NL, Haynes RB, Ramkissoonsingh R; Hedges Team. Developing optimal search strategies for detecting sound clinical prediction studies in MEDLINE. AMIA Annu Symp Proc. 2003;2003:728-32. [PMID: 14728269](https://pubmed.ncbi.nlm.nih.gov/14728269/); [PMCID: PMC1479983](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1479983/).
- Wilczynski NL, Haynes RB; Hedges Team. Developing optimal search strategies for detecting clinically sound causation studies in MEDLINE. AMIA Annu Symp Proc. 2003;2003:719-23. [PMID: 14728267](https://pubmed.ncbi.nlm.nih.gov/14728267/); [PMCID: PMC1480286](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1480286/).
- Wilczynski NL, McKibbon KA, Haynes RB. Enhancing retrieval of best evidence for health care from bibliographic databases: calibration of the hand search of the literature. Stud Health Technol Inform. 2001;84(Pt 1):390-3. [PMID: 11604770](https://pubmed.ncbi.nlm.nih.gov/11604770/).
- Haynes RB, Wilczynski N, McKibbon KA, Walker CJ, Sinclair JC. Developing optimal search strategies for detecting clinically sound studies in MEDLINE. J Am Med Inform Assoc. 1994 Nov-Dec;1(6):447-58. doi: 10.1136/jamia.1994.95153434. [PMID: 7850570](https://pubmed.ncbi.nlm.nih.gov/7850570/); [PMCID: PMC116228](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC116228/).
 #### Medical genetics search filters
 The [medical genetics](#medical-genetics) searches were developed in conjunction with the staff of GeneReviews: Genetic Disease Online Reviews at GeneTests, University of Washington, Seattle.
 Category PubMed equivalent Diagnosis (Diagnosis AND Genetics) Differential Diagnosis (Differential Diagnosis[MeSH] OR Differential Diagnosis[Text Word] AND Genetics) Clinical Description (Natural History OR Mortality OR Phenotype OR Prevalence OR Penetrance AND Genetics) Management (therapy[Subheading] OR treatment[Text Word] OR treatment outcome OR investigational therapies AND Genetics) Genetic Counseling (Genetic Counseling OR Inheritance pattern AND genetics) Molecular Genetics (Medical Genetics OR genotype OR genetics[Subheading] AND genetics) Genetic Testing (DNA Mutational Analysis OR Laboratory techniques and procedures OR Genetic Markers OR diagnosis OR testing OR test OR screening OR mutagenicity tests OR genetic techniques OR molecular diagnostic techniques AND genetics) Medical Genetics ((Diagnosis AND genetics) OR (Differential Diagnosis[MeSH] OR Differential Diagnosis[Text Word] AND genetics) OR (Natural History OR Mortality OR Phenotype OR Prevalence OR Penetrance AND genetics) OR (therapy[Subheading] OR treatment[Text Word] OR treatment outcome OR investigational therapies AND genetics) OR (Genetic Counseling OR Inheritance pattern AND genetics) OR (Medical Genetics OR genotype OR genetics[Subheading] AND genetics) OR (DNA Mutational Analysis OR Laboratory techniques and procedures OR Genetic Markers OR diagnosis OR testing OR test OR screening OR mutagenicity tests OR genetic techniques OR molecular diagnostic techniques AND genetics)) ### Computation of similar articles
 The neighbors of a document are those documents in the database that are the most similar to it. The similarity between documents is measured by the words they have in common, with some adjustment for document lengths. To carry out such a program, one must first define what a word is. For us, a word is basically an unbroken string of letters and numerals with at least one letter of the alphabet in it. Words end at hyphens, spaces, new lines, and punctuation. The 132 common, but uninformative, words (also known as stopwords) are eliminated from processing at this stage. Next, a limited amount of stemming of words is done, but no thesaurus is used in processing. Words from the abstract of a document are classified as text words. Words from titles are also classified as text words, but words from titles are added in a second time to give them a small advantage in the local weighting scheme. MeSH terms are placed in a third category, and a MeSH term with a subheading qualifier is entered twice, once without the qualifier and once with it. If a MeSH term is starred (indicating a major concept in a document), the star is ignored. These three categories of words (or phrases in the case of MeSH) comprise the representation of a document. No other fields, such as Author or Journal, enter into the calculations.
 Having obtained the set of terms that represent each document, the next step is to recognize that not all words are of equal value. Each time a word is used, it is assigned a numerical weight. This numerical weight is based on information that the computer can obtain by automatic processing. Automatic processing is important because the number of different terms that have to be assigned weights is close to two million for this system. The weight or value of a term is dependent on three types of information: 1) the number of different documents in the database that contain the term; 2) the number of times the term occurs in a particular document; and 3) the number of term occurrences in the document. The first of these pieces of information is used to produce a number called the global weight of the term. The global weight is used in weighting the term throughout the database. The second and third pieces of information pertain only to a particular document and are used to produce a number called the local weight of the term in that specific document. When a word occurs in two documents, its weight is computed as the product of the global weight times the two local weights (one pertaining to each of the documents).
 The global weight of a term is greater for the less frequent terms. This is reasonable because the presence of a term that occurred in most of the documents would really tell one very little about a document. On the other hand, a term that occurred in only 100 documents of one million would be very helpful in limiting the set of documents of interest. A word that occurred in only 10 documents is likely to be even more informative and will receive an even higher weight.
 The local weight of a term is the measure of its importance in a particular document. Generally, the more frequent a term is within a document, the more important it is in representing the content of that document. However, this relationship is saturating, i.e., as the frequency continues to go up, the importance of the word increases less rapidly and finally comes to a finite limit. In addition, we do not want a longer document to be considered more important just because it is longer; therefore, a length correction is applied. This local weight computation is based on the Poisson distribution and the formula can be found in [Lin J and Wilbur WJ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2212667/).
 The similarity between two documents is computed by adding up the weights (local wt1 × local wt2 × global wt) of all of the terms the two documents have in common. This provides an indication of how related two documents are. The resultant score is an example of a vector score. Vector scoring was originated by Gerard Salton and has a long history in text retrieval. The interested reader is referred to Salton, Automatic Text Processing, Reading, MA: Addison-Wesley, 1989 for further information on this topic. Our approach differs from other approaches in the way we calculate the local weights for the individual terms. Once the similarity score of a document in relation to each of the other documents in the database has been computed, that document's neighbors are identified as the most similar (highest scoring) documents found. These closely related documents are pre-computed for each document in PubMed so that when you select Similar articles, the system has only to retrieve this list. This enables a fast response time for such queries.
 ### Journal lists
 PubMed journals
 - [Uncompressed](https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt)
- [GNU zip](https://ftp.ncbi.nih.gov/pubmed/J_Medline.gz)
- [UNIX Compress](https://ftp.ncbi.nih.gov/pubmed/J_Medline.Z)
- [PKZIP](https://ftp.ncbi.nih.gov/pubmed/J_Medline.zip)
 PubMed and NCBI molecular biology database journals
 - [Uncompressed](https://ftp.ncbi.nih.gov/pubmed/J_Entrez.txt)
- [GNU zip](https://ftp.ncbi.nih.gov/pubmed/J_Entrez.gz)
- [UNIX Compress](https://ftp.ncbi.nih.gov/pubmed/J_Entrez.Z)
- [PKZIP](https://ftp.ncbi.nih.gov/pubmed/J_Entrez.zip)
 Back to Top

