#!/bin/bash
. /etc/profile
#myzsh=${ZSH_NAME:+zsh};
#mybash=${BASH_VERSION:+bash};
#myksh=${PS3:+ksh};
#mysh=${myzsh:-${mybash:-${myksh:-sh}}};
#unset myzsh mybash myksh;
#module load R-Project/3.1.2;
module available;
export PATH=/afs/cad/linux/R-3.1.2/bin:$PATH
Rscript PAEAtest.R;
