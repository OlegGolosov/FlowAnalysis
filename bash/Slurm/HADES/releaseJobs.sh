squeue -u ogolosov | grep JobHeld | awk '{print $1}' | xargs -l1 scontrol release

