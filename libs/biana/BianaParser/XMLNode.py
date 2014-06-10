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


class XMLNode(object):

    def __init__(self, name, attrs):

        self.name = name
        self.attrs = attrs
        self.childs = {}
        self.value = []


    def addChild(self, XMLNodeChild):
        self.childs.setdefault(XMLNodeChild.name.lower(), []).append(XMLNodeChild)

    def getValue(self):
        return "".join(self.value)

    def addValue(self, value):
        self.value.append(value)

    def getChilds(self, child_name=None):
        #print [ y.name for x in self.childs.values() for y in x ]
        if child_name is None:
            c = []
            [ c.extend(x) for x in self.childs.values() ]
            return c
        return self.childs[child_name.lower()]

    def getChild(self, child_name):
        if( len(self.childs[child_name.lower()])!=1 ):
            #raise ValueError("You are trying to obtain a single child when it has %s childs" %(len(self.childs[child_name.lower()]) ) )
            return None
        return self.childs[child_name.lower()][0]

