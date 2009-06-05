

######################################################   HTML

PREFIX ?= /usr/local
WEB_DESTDIR ?= /tmp/localweb

ASCIIDOC ?= asciidoc --unsafe --backend=xhtml11 --attribute icons --attribute iconsdir=/images/icons

CXXFLAGS += -std=c++0x -I/usr/local/include -I$(PREFIX)/include/lopti

show: localweb
	firefox $(WEB_DESTDIR)/index.html

index.html: index.txt
	$(ASCIIDOC)  -o $@  $<

localweb: 	 index.html *.gif example*.cc
	mkdir -p  $(WEB_DESTDIR)
	cp -uv $^ $(WEB_DESTDIR)

clean:
	rm -f *.html /tmp/localweb
	git clean -df

######################################################   EXAMPLES

example:	

example-gsl:	LDFLAGS=-lgsl

example-newuoa:	LDFLAGS=-llopti -lgfortran

