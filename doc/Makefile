

PREFIX ?= /usr/local
#NAME = lopti
WEB_DESTDIR ?= /tmp/localweb
#WEB_DESTDIR := $(WEB_DESTDIR)/$(NAME)

ASCIIDOC ?= asciidoc --unsafe --backend=xhtml11 --attribute icons --attribute iconsdir=/images/icons

CXXFLAGS += -std=c++0x -I/usr/local/include -I$(PREFIX)/include/lopti

show: localweb
	firefox $(WEB_DESTDIR)/index.html

example:

example-gsl:	LDFLAGS=-lgsl
example-newuoa:	LDFLAGS=-llopti

index.html: index.txt
	$(ASCIIDOC)  -o $@  $<

localweb: 	 index.html *.gif example*.cc
	mkdir -p  $(WEB_DESTDIR)
	cp -uv $^ $(WEB_DESTDIR)

clean:
	rm -f *.html
	git clean -df
