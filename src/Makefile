# FRAME3DD: 
# Static and dynamic structural analysis of 2D & 3D frames and trusses
# with elastic and geometric stiffness.
# ---------------------------------------------------------------------------
# http://frame3dd.sourceforge.net/
# ---------------------------------------------------------------------------
# Copyright (C) 1992-2009  Henri P. Gavin
#
#    FRAME3DD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FRAME3DD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
##  
#  @file
#            a simple Makefile for Frame3DD
## 

SHELL = /bin/sh

CC = gcc
# CFLAGS = -Wall -c -O 
  CFLAGS = -c -O 
LDFLAGS = -lm
SOURCES =  main.c frame3dd.c frame3dd_io.c coordtrans.c eig.c HPGmatrix.c HPGutil.c NRutil.c 
OBJECTS = $(SOURCES:.c=.o)
EXECUTEABLE = frame3dd

all : $(SOURCES) $(EXECUTEABLE)

$(EXECUTEABLE) : $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@
  
.c.o : 
	$(CC) $(CFLAGS) $< -o $@

clean :
	rm -rf *.o

install :
	mv frame3dd /usr/local/bin/.


