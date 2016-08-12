#!/bin/sh
SESSION="hotspotr"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n R
tmux send-keys -t $SESSION:1 'vim README.Rmd' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/neutral-hotspots.R' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/neutral-hotspots-ntests.R' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/fit-hotspot-model.R' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/rs-dist-diff.R' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/p-values.R' C-m
#tmux send-keys -t $SESSION:1 '2gt'

tmux new-window -t $SESSION:2 -k -n cpp
tmux send-keys -t $SESSION:2 'cd ./src' C-m
tmux send-keys -t $SESSION:2 'vim ives.cpp' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe neutral-hotspots.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe neutral-hotspots.cpp' C-m
tmux split-window -h
tmux send-keys -t $SESSION:2 'cd ./src' C-m
tmux send-keys -t $SESSION:2 'vim ac-stats.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe ac-stats.cpp' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:3 -n makefile
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim tmux-start.bash' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe makefile' C-m
tmux split-window -h
tmux send-keys -t $SESSION:3 'xdg-open README.html &' C-m
tmux send-keys -t $SESSION:3 'git status -uno' C-m
tmux select-pane -t 0

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
