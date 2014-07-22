#!/bin/bash
#
# Top level Makefile for MAPS                   DAM Sep 27 2007
#

BASEDIR = $(PWD)
SUBDIRS = novas-c201 utilities maps_im2uv visgen maps2uvfits
BINEXECUTABLES = maps2uvfits maps_im2uv visgen

install:
	@for subdir in $(SUBDIRS); do \
	  echo "Installing all in $$subdir"; \
	  cd $$subdir && $(MAKE) install && cd $(BASEDIR); \
	done

clean:
	@for subdir in $(SUBDIRS); do \
	  echo "Cleaning object files in $$subdir"; \
	  cd $$subdir && $(MAKE) clean && cd $(BASEDIR); \
	done

purge:
	@echo 'Cleaning remporary emacs *~ and vi #*# files'; \
	rm -f `find ./ -name '*~'` && rm -f `find ./ -name '#*#'`; \
	echo 'Cleaning executables in /bin'; \
	cd $(BASEDIR)/bin && rm -f $(BINEXECUTABLES)  && cd $(BASEDIR); \
	echo 'Cleaning libraries in lib'; \
	cd $(BASEDIR)/lib && rm -f *.a  && cd $(BASEDIR); \
