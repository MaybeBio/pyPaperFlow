# PubMed platform demo: search → fetch metadata → metadata display.
#
# Render from the repository root (needs vhs + ttyd + ffmpeg + a Chromium):
#     http_proxy= https_proxy= vhs figs/pubmed.tape
#
# `paperflow` must be on PATH inside the recorded shell. The demo runs in a
# scratch directory so nothing is written into the repository.
#
# NCBI requires a contact email (--email). This tape uses a placeholder so no
# real address appears on screen; substitute your own when running it yourself.
# Adding --api-key raises the NCBI rate limit from 3 to 10 requests/second.
#
# The query is deliberately narrow: pubmed-search retrieves *every* matching
# PMID (the -n option is the batch size, not a result cap), and pubmed-meta
# fetches linked data per paper, so a broad query takes minutes.
Output "figs/pubmed.gif"

Set Shell "bash"
Set FontFamily "DejaVu Sans Mono"
Set FontSize 16
Set Width 1000
Set Height 620
Set Padding 16
Set Framerate 15
Set TypingSpeed 40ms
Set PlaybackSpeed 2.6
Set Theme "Dracula"
Set CursorBlink false

# Hidden setup: scratch dir + plain prompt + clean screen.
Hide
Type "mkdir -p /tmp/paperflow-pubmed && cd /tmp/paperflow-pubmed && rm -rf out papers && PS1='$ '"
Enter
Type "clear"
Enter
Show
Sleep 400ms

Type "paperflow pubmed-search 'AlphaFold3 AND covalent AND inhibitor AND 2025[dp]' --email you@example.com -o ./out"
Enter
Sleep 8s

Sleep 400ms
Type "paperflow pubmed-meta --file ./out/pubmed_searched_ids_*.txt --email you@example.com -s ./papers -o ./out"
Enter
Sleep 24s

Sleep 400ms
Type "jq -r '.identity.title' ./out/pubmed/*/*/*_meta.json"
Enter
Sleep 3s

Sleep 400ms
Type "head -n 14 $(ls ./out/pubmed/*/*/*_meta.json | head -1)"
Enter
Sleep 3s
