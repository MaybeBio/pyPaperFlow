# ChemRxiv platform demo: search → fetch metadata → metadata display.
#
# Render from the repository root (needs vhs + ttyd + ffmpeg + a Chromium):
#     http_proxy= https_proxy= vhs figs/chemrxiv.tape
#
# `paperflow` must be on PATH inside the recorded shell. The demo runs in a
# scratch directory so nothing is written into the repository.
Output "figs/chemrxiv.gif"

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
Type "mkdir -p /tmp/paperflow-chem && cd /tmp/paperflow-chem && rm -rf out && PS1='$ '"
Enter
Type "clear"
Enter
Show
Sleep 400ms

Type "paperflow chemrxiv-search 'machine learning potential' -n 5 -o ./out"
Enter
Sleep 12s

Sleep 400ms
Type "paperflow chemrxiv-fetch --file ./out/searched_chemrxiv_ids.txt -o ./out"
Enter
Sleep 30s

Sleep 400ms
Type "jq -r .title ./out/chemrxiv/*/*/*.json"
Enter
Sleep 3s

Sleep 400ms
Type "head -n 9 $(ls ./out/chemrxiv/*/*/*.json | head -1)"
Enter
Sleep 3s
