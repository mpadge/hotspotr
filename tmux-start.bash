#!/bin/sh
SESSION="hotspotr"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim README.Rmd' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/test2d.R' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe src/neutral2d.cpp' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe R/run_tests.R' C-m
#tmux send-keys -t $SESSION:1 '2gt'
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -n makefile
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim tmux-start.bash' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe makefile' C-m
tmux split-window -h
tmux send-keys -t $SESSION:2 'git status -uno' C-m
tmux select-pane -t 0

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
