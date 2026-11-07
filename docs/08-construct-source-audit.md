# Exchange-construct source audit

Reviewed 2026-11-07 in a bounded public-source search. **Exact complete-genome mapping remains unverified for both whole-gp17 exchanges.** Four phenotype observations therefore remain abstentions. This is not evidence that a public sequence does not exist.

The [original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4785837/) distinguishes whole-gp17 and C-terminal exchange constructs and reports checks of the heterologous region. Those statements do not establish a public complete-genome record for each experimental stock.

The [PubMed-to-nuccore lookup](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=nuccore&id=26973885&retmode=json) returned no linked nucleotide records. Narrow title and construct-name searches did not establish a mapping; deposited names may differ. The related records [LC730324.1, T3_AGU](https://www.ncbi.nlm.nih.gov/nuccore/LC730324.1) and [LC730325.1, T7_AGU](https://www.ncbi.nlm.nih.gov/nuccore/LC730325.1) describe 2022 material and have not been linked to the 2015 exchanges. They cannot be substituted merely because their names contain T3 or T7.

## Access limits

The paper's linked PDF and XLSX supplements returned reCAPTCHA HTML rather than the actual attachments in this audit. The [Europe PMC supplement endpoint](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4785837/supplementaryFiles) reported that the article was not open access, and the publisher returned 403. Supplement contents were not inspected. The [machine-readable record](../software/tipcraft/benchmarks/retro-t7-t3/construct-source-audit.json) lists URLs, findings and limits for review.

## Consequence for evaluation

Preserve the four abstentions until a complete-genome record is explicitly mapped to each exact experimental construct. Native reference genomes, other recombinant isolates, short patent fragments and inferred reconstructions do not establish that correspondence. No new sequences were built, and no authors or classmates were contacted during this audit.

The report now accepts `--manifest` to verify the frozen input hash, parameters, candidate sets, exclusion declarations and both adjacent upstream likelihood matrices. Excluded constructs with scores and runs declaring label access are rejected. Editing JSON scores alone cannot create a verified covered case. These checks establish artifact consistency; they do not independently authenticate a deposited sequence or experimental stock.
