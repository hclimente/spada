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

import traceback
import socket
import sys

class OutBianaInterface(object):

    outmethod = None
    out_format = "xml"

    def connect_to_socket(port,server="127.0.0.1"):
        """
        Changes the default outmethod to a socket
        It should not be used by users, only graphical interface should use this method
        """
        
        serverHost = server
        serverPort = port
        
        try:
            # Start socket communication
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)    # create a TCP socket
            s.connect((serverHost, serverPort)) # connect to server on the port
            s.send("<?xml version=\"1.0\"?>\r\n")
            s.send("<biana_to_gui>")
            OutBianaInterface.outmethod = s.send
        except:
            print "Impossible to connect to port %s" %port
            traceback.print_exc()
            pass

    connect_to_socket = staticmethod(connect_to_socket)

    def set_outmethod(method):
        OutBianaInterface.outmethod = method

    set_outmethod = staticmethod(set_outmethod)

    def send_data(message):
        if OutBianaInterface.outmethod is not None:
            OutBianaInterface.outmethod(message)
            #sys.stderr.write("M: %s\n"%message)
        #sys.stderr.write("M: %s\n"%message)

    send_data = staticmethod(send_data)

    def send_process_message(message):
        if OutBianaInterface.out_format == "xml":
            OutBianaInterface.send_data("<start_biana_process process=\"%s\" />" %message)
        else:
            OutBianaInterface.send_data("%s\n" %message)

    send_process_message = staticmethod(send_process_message)

    def send_end_process_message():
        if OutBianaInterface.out_format == "xml":
            OutBianaInterface.send_data("<end_biana_process />")

    send_end_process_message = staticmethod(send_end_process_message)

    def send_error_notification(message,error):
        if OutBianaInterface.out_format == "xml":
            # To avoid parser error in case error line contains '<', '>' or '&'
            error = "<![CDATA[%s]]>" % error
            OutBianaInterface.send_data("<error_notification value=\"%s\">%s</error_notification>" %(message,error))
            #sys.stderr.write("M: %s\nE: %s\n" % (message, error))
        else:
            OutBianaInterface.send_data("\n!%s:\n\n%s\n" %(message, error))

        if OutBianaInterface.outmethod is None:
            sys.stderr.write(message)
            sys.stderr.write("\n")
            sys.stderr.write(error)
            sys.stderr.write("\n")

    send_error_notification = staticmethod(send_error_notification)

    def send_info_message(message):
        if OutBianaInterface.out_format == "xml":
            message = "<![CDATA[%s]]>" % message
            OutBianaInterface.send_data("<info_message>%s</info_message>" %(message))
        else:
            OutBianaInterface.send_data("%s" %message)

    send_info_message = staticmethod(send_info_message)


    def close():
        OutBianaInterface.send_data("</biana_to_gui>")


    close = staticmethod(close)
