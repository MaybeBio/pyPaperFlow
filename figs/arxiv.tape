# arXiv platform demo: search → fetch metadata + PDFs → metadata display.
#
# Render from the repository root (needs vhs + ttyd + ffmpeg + a Chromium):
#     http_proxy= https_proxy= vhs figs/arxiv.tape
#
# `paperflow` must be on PATH inside the recorded shell. The demo runs in a
# scratch directory so nothing is written into the repository.
Output "figs/arxiv.gif"

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
Type "mkdir -p /tmp/paperflow-arxiv && cd /tmp/paperflow-arxiv && rm -rf out papers && PS1='$ '"
Enter
Type "clear"
Enter
Show
Sleep 400ms

Type "paperflow arxiv-search 'protein structure prediction' -n 5 -o ./out"
Enter
Sleep 14s

Sleep 400ms
Type "paperflow arxiv-fetch --file ./out/searched_arxiv_ids.txt -s ./papers -o ./out"
Enter
Sleep 35s

Sleep 400ms
Type "jq -r .title ./out/arxiv/*/*/*.json"
Enter
Sleep 3s

Sleep 400ms
Type "find ./out/arxiv -name '*.pdf' | head -5"
Enter
Sleep 3s

Sleep 400ms
Type "head -n 9 $(ls ./out/arxiv/*/*/*.json | head -1)"
Enter
Sleep 3s
