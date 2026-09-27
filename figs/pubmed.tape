# PubMed platform demo: search → fetch metadata → metadata display.
#
# Render from the repository root (needs vhs + ttyd + ffmpeg + a Chromium):
#     http_proxy= https_proxy= vhs figs/pubmed.tape
#
# NCBI requires a contact email (--email). Replace the address below with your
# own before recording — it is visible on screen. Adding --api-key raises the
# rate limit from 3 to 10 requests/second and is recommended.
#
# NOT YET RECORDED: NCBI's E-utilities backend is currently returning HTTP 500.
# This tape is ready to run once it recovers.
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

Type "paperflow pubmed-search 'AlphaFold protein structure' --email luxunisgod123@gmail.com -n 5 -o ./out"
Enter
Sleep 20s

# pubmed-search names its id file with a timestamp, so glob the newest one.
Sleep 400ms
Type "paperflow pubmed-meta --file $(ls -t ./out/pubmed_searched_ids_*.txt | head -1) --email luxunisgod123@gmail.com -s ./papers -o ./out"
Enter
Sleep 40s

Sleep 400ms
Type "jq -r .title ./out/pubmed/*/*/*.json"
Enter
Sleep 3s

Sleep 400ms
Type "head -n 12 $(ls ./out/pubmed/*/*/*.json | head -1)"
Enter
Sleep 3s
