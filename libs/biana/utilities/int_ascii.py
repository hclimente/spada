"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#PRIMARY METHODS

# For signed and unsigned values

def int_to_ascii(max_bytes,value):
    value = int(value)
    if value<0:
        value = -1*value | pow(2,max_bytes*8-1)
    temp = [ '\x00' for x in xrange(max_bytes) ]
    del x
    pos = max_bytes-1
    while value>255:
        temp[pos] = chr(value%255)
        pos-=1
        value/=255
    if value==255: temp[pos] = chr(255)
    else: temp[pos] = chr(value%255)
    return "".join(temp)

def ascii_to_int(ascii):
    value = 0
    pos = len(ascii)-1
    for x in xrange(len(ascii)):
        value += ord(ascii[x])*pow(255,pos)
        pos-=1
    if ord(ascii[0])&128:
        return -1* value ^ pow(2,len(ascii)*8-1)
    else:
        return value


# For unsigned values

def unsigned_int_to_ascii(max_bytes,value):

    value = int(value)
    temp = [ '\x00' for x in xrange(max_bytes) ]
    del x
    pos = max_bytes-1
    while value>255:
        temp[pos] = chr(value%255)
        pos-=1
        value/=255
    if value==255: temp[pos] = chr(255)
    else: temp[pos] = chr(value%255)

    return "".join(temp)


def ascii_to_unsigned_int(ascii):

    value = 0

    pos = len(ascii)-1
    
    for x in xrange(len(ascii)):
        value += ord(ascii[x])*pow(255,pos)
        pos-=1

    return value



# DERIVED METHODS

def unsigned_float_to_ascii(bytes,dec_positions,value):

    return unsigned_int_to_ascii(bytes,float(value)*pow(10,dec_positions))
    
def ascii_to_unsigned_float(dec_positions, ascii):

    return ascii_to_unsigned_int(ascii)/float(pow(10,dec_positions))

def float_to_ascii(bytes,dec_positions,value):

    return int_to_ascii(bytes,float(value)*pow(10,dec_positions))

def ascii_to_float(dec_positions,ascii):

    return ascii_to_int(ascii)/float(pow(10,dec_positions))



# ITERATIVE METHODS


def unsigned_int_list_to_ascii(bytes,list):

    return "".join( [unsigned_int_to_ascii(bytes,x) for x in list] )

def ascii_to_unsigned_int_list(bytes,ascii):

    return [ ascii_to_unsigned_int(ascii[x:x+bytes]) for x in xrange(0,len(ascii),bytes)]

def unsigned_float_list_to_ascii(int_bytes,dec_bytes,list):

    return "".join( [unsigned_float_to_ascii(int_bytes,dec_bytes,x) for x in list] )

def ascii_to_unsigned_float_list(bytes,dec_bytes,ascii):

    return [ ascii_to_unsigned_float(dec_bytes,ascii[x:x+bytes]) for x in xrange(0,len(ascii),bytes)]


def int_list_to_ascii(bytes,list):

    return "".join( [int_to_ascii(bytes,x) for x in list] )

def ascii_to_int_list(bytes,ascii):

    return [ ascii_to_int(ascii[x:x+bytes]) for x in xrange(0,len(ascii),bytes)]

def float_list_to_ascii(int_bytes,dec_bytes,list):

    return "".join( [float_to_ascii(int_bytes,dec_bytes,x) for x in list] )

def ascii_to_float_list(bytes,dec_bytes,ascii):

    return [ ascii_to_float(dec_bytes,ascii[x:x+bytes]) for x in xrange(0,len(ascii),bytes)]
