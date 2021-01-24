# To update, edit /opt/local/etc/skel

# uncomment the following to activate bash-completion:
[ -f /etc/profile.d/bash-completion ] && source /etc/profile.d/bash-completion

###############################################################################
#
# Specific module account settings:
#
###############################################################################
MODULE_INIT=/opt/Modules/init/bash
if [ -f $MODULE_INIT ]; then
	source /opt/Modules/init/bash
	module load ics-default
else
	function module {
		if [ -z $1 ]; then
			echo "Modules are not available on this system"
		fi
	}
fi

###############################################################################
#
# Variables
#
###############################################################################
# Test for an interactive shell.  There is no need to set anything
# past this point for scp and rcp, and it's important to refrain from
# outputting anything in those cases.
if [[ $- != *i* ]]; then
        # Shell is non-interactive.  Be done now
	return
fi


#add the printer name!
#export PRINTER=name
#export LPRINTER=name

# history files
# change the size 
export HISTFILE="~/.bash_history"
export HISTSIZE=1000
export HISTFILESIZE=10000

# Change the window itle of X terminals 
PROMPT_COMMAND='echo -ne "\033]0;${USER}@${HOSTNAME%%.*}:${PWD/$HOME/~}\007"'

#prompt
BLACK="\[\033[0;30m\]"
BLUE="\[\033[0;34m\]"
GREEN="\[\033[0;32m\]"
GREY="\[\033[01;30m\]"
CYAN="\[\033[0;36m\]"
RED="\[\033[01;31m\]"
PURPLE="\[\033[0;35m\]"
BROWN="\[\033[0;33m\]"
LGREY="\[\033[0;37m\]"
Yellow="\[\033[0;33m\]"
White="\[\033[0;37m\]"
NORMAL="\[\033[00m\]"

#PS1="${GREY}\u@\h ${GREEN}\t ${RED}\w \n\$ ${NORMAL}"
#comment this PS1 and use above if you want color
PS1="\u@\h \t \w \n\$ "

# The default for PS2 is > which may be mistaken for a re-direct
export PS2="\\"

#General Settings
set noclobber
# Let me know when a background job has finished the moment it finishes.
set notify

# limit the core size to 1MB (2000 512-byte blocks)
# ulimit -c 2000

# Make sure files are NOT world readable
umask 077

THEOS=`uname`
THEREV=`uname -r`
RUID=`/usr/[ub][ci][bn]/whoami`

###############################################################################
#
# Specific alias account settings:
# These aliases are for all operating systems
#
###############################################################################

#navigation
alias up="cd .."

#graphical
alias xterm="xterm -rv -sb &"

# 
alias s="suspend"
alias 1="fg %1"
alias 2="fg %2"

#ph alias'
alias pha="ph alias=\!* return all"
alias phn="ph name=\!* return all"
alias phe="ph ext=\!* return all"
alias phi="ph id=\!*"

###############################################################################
#
# Specific Linux settings:
#
###############################################################################

if test "$THEOS" == "Linux" ; then
	#settings
	# colors for ls, etc. 
	if test -f "/etc/DIR_COLORS"; then
		eval `dircolors -b /etc/DIR_COLORS`
	fi
	#aliases
	alias konsole="konsole --noxft"
	alias d="ls --color"
	alias cp="cp -iv"
	alias rm="rm -iv"
	alias mv="mv -iv"
	alias ls="ls -F --color=auto --human-readable --almost-all"
	alias ll="ls -l -F --color=auto --human-readable --almost-all"
	alias df="df -h"

###############################################################################
#
# Specific Solaris settings:
#
###############################################################################

elif test "$THEOS" == "SunOS" ; then
	#settings
	TERM=vt100
	#aliases
	alias ls="ls -F"
	alias ll="ls -la"
	alias cp="cp -i"
	alias mv="mv -i"
	alias rm="rm -i"
	alias reboot="echo Sure?"
fi

