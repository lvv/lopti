HOMEDIR ?= /home/lvv/p/volnitsky.com/
include $(HOMEDIR)/include.mk

######################################################   EXAMPLES

PREFIX ?= /usr/local
CXXFLAGS += -std=c++0x -I/usr/local/include -I$(PREFIX)/include/lopti

example:	

example-gsl:	LDFLAGS=-lgsl

example-newuoa:	LDFLAGS=-llopti -lgfortran

