# medRxiv platform demo: search → fetch metadata → metadata display.
#
# Render from the repository root (needs vhs + ttyd + ffmpeg + a Chromium):
#     http_proxy= https_proxy= vhs figs/medrxiv.tape
#
Output "figs/medrxiv.gif"

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
Type "mkdir -p /tmp/paperflow-med && cd /tmp/paperflow-med && rm -rf out && PS1='$ '"
Enter
Type "clear"
Enter
Show
Sleep 400ms

Type "paperflow medrxiv-search 'sleep disordered breathing' -n 5 -o ./out"
Enter
Sleep 14s

Sleep 400ms
Type "paperflow medrxiv-fetch --file ./out/searched_medrxiv_ids.txt -o ./out"
Enter
Sleep 30s

Sleep 400ms
Type "jq -r .title ./out/medrxiv/*/*/*.json"
Enter
Sleep 3s

Sleep 400ms
Type "head -n 9 $(ls ./out/medrxiv/*/*/*.json | head -1)"
Enter
Sleep 3s
