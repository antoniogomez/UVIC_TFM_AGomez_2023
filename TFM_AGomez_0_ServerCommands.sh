# YAMABUKI SERVER
## Interactive Session

### after typing password to enter the yamabuki, then: 
  
module load tools/tmux/2.3 # "tmux" works like "screen" (if it's the first time then type module av tools/tmux/2.3, close the session and open it again directly typing the next line)

tmux new -s InteractiveJob # or any other name you want to give to the interactive session

srun --job-name "InteractiveJob" --partition=short --mem-per-cpu 8GB --pty bash # if "short" then 12 hours of interactive session, if "long" then 72 hours; you can also change the cpu needed


``
## In Yamabuki server
###  Interactive Session
interactive job of maximum 2h (connection is slow)
```{bash}
# after typing password to enter the yamabuki, then: 

module load tools/tmux/2.3 
# "tmux" works like "screen" (if it's the first time then type module av tools/tmux/2.3, close the session and open it again directly typing the next line)

tmux new -s InteractiveJob 
# or any other name you want to give to the interactive session
```

- option 1: this command (cpus from 2 to 8)
```{bash}
srun --job-name "InteractiveJob" --cpus-per-task 2 --mem-per-cpu 4096  --time 2:00:00 --pty bash
```
o
```{bash}
srun --job-name "InteractiveJob" --partition=short --mem-per-cpu 8GB --pty bash
# if "short" then 12 hours of interactive session, if "long" then 72 hours; you can also change the cpu needed
```

- option 2: this command with .sh
```{bash}
cd /PROJECTES/BISC_OMICS
sh interactive.sh
```

Look at modules available 
```{bash}
module av
```

look for a module
```{bash}
module spider <nom del que vols buscar>
  ```

load plink2 module
```{bash}
module load bio/PLINK/1.9b_6.21-x86_64 #There is another versions: 2.00a3.6-GCC-11.3.0 and 2.00a2.3_x86_6
#module load bio/PLINK/2.00a3.6-GCC-11.3.0
```

### R
```{bash}
# these three lines are used to start working with R
module load system/OpenSSL/1.1
module load lang/R/4.1.2-foss-2021b
R   

# to close the R session + tmux session + interactive session: 
quit()
save workspace y/n/ ? n
tmux kill-session -t InteractiveJob
exit

# if you don't want to close the session, then just close mobaxterm and 
# to connect again to the same session (within the 12 or X hours requested) then type: 
tmux ls #(shows all the sessions)
tmux a -t InteractiveJob #(reattaches the session indicated)

# to perform a dettach and go back to the "sessió mare", then type "ctrl+b+d". 
```
