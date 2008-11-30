
#include ../Makefile

PREFIX ?= /usr/local
DESTDIR ?= /tmp/html-lopti
ASCIIDOC ?= asciidoc

CXXFLAGS += -I$(PREFIX)/include/lopti

show: install
	firefox $(DESTDIR)/index.html

example:

example-gsl:	LDFLAGS=-lgsl
example-newuoa:	LDFLAGS=-llopti

index.html: index.txt
	$(ASCIIDOC)  -o $@  $<

install: index.html *.gif example*.cc
	mkdir -p  $(DESTDIR)
	cp -uv $^ $(DESTDIR)

clean:
	rm -f *.html
	git clean
