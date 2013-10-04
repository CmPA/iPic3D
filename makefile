# Convenience makefile to call scripts

help:
	scripts/ipic help

tags: retags

retags:
	scripts/ipic ctags

#monitor:
#	less +F data/ConservedQuantities.txt
