##############################################################################
##                   Copyright (C) 1997 CEA
##############################################################################
##    UNIT
##
##    Version: 1.10
##
##    Author: Jean-Luc Starck
##
##    Date:  98/03/26
##    
##    File:  Makefile
## 
##############################################################################
##
## MAKEFILE
## --------
##  Le but de ce $(MAKE)file est de gerer les librairies du package
##  il permet de compiler les differentes librairies, de compter
##  le nomber de lignes de programmes, de voir les fichiers pris
##  en edition, d'imprimer les fichier ...
##  les differentes options sont decrites ci-dessous       
##
## OPTION du Makefile
## -------------------
##
## Compilation de tout le package + Outils
##
##    . $(MAKE) $(MAKE)files 
##           lance la compilation de toutes les librairies
##           
##    . $(MAKE) listing 
##           imprime tous les sources en utilisant le programme a2ps
##
##    . $(MAKE) sccsget 
##           extrait tous les sources des directories SCCS
##
##    . $(MAKE) sccsinfo 
##           indique quels sont les fichiers pris en edition
##
##    . $(MAKE) nb_lines
##           compte le nombre de ligne de programmes qui ne sont
##           pas des lignes blanches
##
##    . $(MAKE) nb_nc_lines
##           compte le nombre de ligne de programmes qui ne sont
##           pas des lignes blanches et qui ne sont pas des commentaires
##
##    . $(MAKE) semicols
##           compte le nombre de ligne de programmes qui 
##           finissent pas un point virgule (instruction)
##    . $(MAKE) lines
##           appelle $(MAKE) nb_lines, nb_nc_lines et semicols
## 
##    . $(MAKE) wc
##           compte le nombre de ligne (lignes blanches + commentaire
##                                       + instruction)
##
##    . clean
##           netoie les librairies
##              suppression des objects modules + executables + sccs clean
##
##    make EXP=/Users/starck/Main/ISAP/idl/TTV2  export    
##         Export toute la librairie IDL
## 
##   make  EXP=/Users/starck/Main/ISAP/idl/TTV2   sccsexport
##         Export toute la librairie IDL, en faisant les liens SCCS (i.e. export pour developpeur).
##
##   make  VERS=2.3 sccsedit
##
##############################################################################

PROJMRS=$(ISAP)
PROJ=$(PROJMRS)/idl

GET = sccs get
LIB =  MISC  MRS  MRSP  Obsolete  COSMO/CMB COSMO/Galaxies COSMO/ISW COSMO/Weaklensing Poisson  STAT  extern/MR1  TEST
EE=$(EXP)
VE=$(VERS)
EXTLIB= MASTER astron icosmo_v1.1 mapview piolib

# BIN = $(PROJROOT)/src
# OWN_CODE = $(LIBRARIES) $(BIN)
# 	    echo === $$i to $(DD) ===;\
# 	    cp -r $$i/*.pro $(DD);\

all: sccsget sccsinfo
current: all 
install: export $(MAKE)files


createtree:
	$ echo "== Export to $(EE) === "
	@for i in $(LIB);\
	do (cd $(PROJ); \
 	    echo === $$i ===;\
	    mkdir -p $(EE)/idl/$$i;\
	    echo);\
	done

copylib:
	@for i in $(EXTLIB);\
	do (cd $(PROJ); \
 	    echo === $$i ===;\
 	    cp -r -f $(PROJ)/extern/$$i  $(EE)/idl/extern;\
 	    echo);\
	done

export: sccsget createtree copylib
	$ echo "== Export to $(EE) === "
	cd $(EE)/idl
	cp isap.pro $(EE)/idl
	cp -r HELP  $(EE)/idl
	@for i in $(LIB);\
	do (cd $(PROJ); \
 	    echo === $$i ===;\
	    echo === $(EE)/idl/$$i ===;\
	    mkdir -p $(EE)/idl/$$i;\
	    cp -f $$i/*.pro $(EE)/idl/$$i;\
	    echo);\
	done

createlink:
	$ echo "== Export to $(EE) === "
	cd $(EE)/idl
	ln -s $(PROJ)/SCCS SCCS
	cp -r HELP  $(EE)/idl
	@for i in $(LIB);\
	do (cd $(EE)/idl/$$i; \
 	    echo === $$i/SCCS to $(EE)/idl/$$i ===;\
	    ln -s $(PROJ)/$$i/SCCS SCCS;\
	    sccs get SCCS;\
	    cd $(EE)/idl/$$i;\
	    echo);\
	done

sccsexport: createtree copylib createlink

sccsinfo:
	@for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo === $$i ===;\
	    sccs info;\
	    echo);\
	done

sccsget:
	@for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo === $$i ===;\
	    sccs get SCCS;\
	    echo);\
	done

sccsedit:
	@for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo === $$i ===;\
	    sccs edit -r$(VE)  SCCS;\
	    echo);\
	done
	
sccsdelget:
	@for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo === $$i ===;\
	    sccs delget 'sccs tell -u';\
	    echo);\
	done

sccsclean:
	@for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo === $$i ===;\
	    sccs clean;\
	    echo);\
	done

nb_lines:
	@ echo "Non blank lines:"
	@ for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo  "$$i: ";\
	    cat *.cc *.h *.java *.pro *.tex 2>bidon | egrep -v '^ *$$|^$$'  | wc -l;\
	    rm -f bidon);\
	done | awk '{t+=$$1;print} END {print "TOTAL: " t}'
	@ echo

nb_nc_lines:
	@ echo "Non blank lines, non comment lines :"
	@ for i in $(OWN_CODE);\
	do (cd $(PROJ)/$$i; \
	    echo "$$i: ";\
	    cat *.pro 2>bidon | egrep -v '^ *$$|^$$|^;*$$|^$$' | wc -l;\
	    rm -f bidon);\
	done | awk '{t+=$$1;print} END {print "TOTAL: " t}'
	@ echo

wc:
	@ echo "Number lines:"
	@ for i in $(LIB);\
	do (cd $(PROJ)/$$i; \
	    echo "$$i: ";\
	    cat *.cc *.h *.java *.pro *.tex 2>bidon | wc -l;\
	    rm -f bidon);\
	done | awk '{t+=$$1;print} END {print "TOTAL: " t}'

lines: nb_nc_lines nb_lines wc



