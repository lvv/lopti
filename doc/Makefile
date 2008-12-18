
#include ../Makefile

PREFIX ?= /usr/local
WEB_DESTDIR ?= /tmp/html-lopti
ASCIIDOC ?= asciidoc

CXXFLAGS += -I$(PREFIX)/include/lopti

show: web_install
	firefox $(WEB_DESTDIR)/index.html

example:

example-gsl:	LDFLAGS=-lgsl
example-newuoa:	LDFLAGS=-llopti

index.html: index.txt
	$(ASCIIDOC)  -o $@  $<

web_install: index.html *.gif example*.cc
	mkdir -p  $(WEB_DESTDIR)
	cp -uv $^ $(WEB_DESTDIR)

clean:
	rm -f *.html
	git clean -df
