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

# Class to manage/handle (create, analyze, operate vs..) networks 
# from userEntities and available interaction data in BianaDB (using BianaDBaccess) 
import sys

class BianaReport(object):
    """
    A class to control BIANA reports
    
    
    """
    
    def __init__(self, title="BIANA report", streamOutputMethod=sys.stdout.write, format="html", url=""):
        """

        title is the title that will be written on top of the report (default is 'BIANA report')
        
        streamOutputMethod is the write method where the report will be written (e.g. file_object.write)
              --> if no streamOutputMethod is given, report is written to stdout

        format is the formatting that will be applied to the report
             - 'html': report in html
             - 'txt': report in raw text

        url is the web address where the interface to BIANA is placed


        """

        self.outmethod = streamOutputMethod

        self.format = format
        
        
        if self.format == 'txt':
            self.outmethod("\n%s\n-----------------------------------------\n" %title)

        elif self.format == "html":
            self.outmethod("<html>\n<head>\n<title>%s</title>\n</head>\n" %(title))
            self.outmethod("<body>")

            # Print the headers of the page: Biana logo, Grib logo, etc etc
            self.outmethod("""<table class="default" border="0" cellpadding="0" cellspacing="0" width="100%">""")
            self.outmethod("""  <tbody>""")
            self.outmethod("""<tr><td>""")
            self.outmethod("""<a href='%s'><img src='%s/images/biana_logo.png' border=0><font size=0></a></td>""" %(url,url))

 
            self.outmethod("""<td class="links" align="center"> """)

            self.outmethod("""<a href="http://sbi.imim.es/staff.html" onmouseover="window.status='Who is who in Structural Bioinformatics Research Laboratory';">People</a> &nbsp;<font style="color: rgb(0, 80, 80);">|</font>&nbsp;""") 

            self.outmethod("""<a href="http://sbi.imim.es/research.html" onmouseover="window.status='Research at our Group';">Research</a> &nbsp;<font style="color: rgb(0, 80, 80);">|</font>&nbsp;""") 

            self.outmethod("""<a href="http://sbi.imim.es/resources.html" onmouseover="window.status='Software developed in our Group';">Resources</a>&nbsp;<font style="color: rgb(0, 80, 80);">|</font>&nbsp;""")

            self.outmethod("""<a href="http://sbi.imim.es/publications.html" onmouseover="window.status='Publications by our Group';">Publications</a> &nbsp;<font style="color: rgb(0, 80, 80);">|</font>&nbsp;""")

            self.outmethod("""<a href="http://sbi.imim.es/links.html" onmouseover="window.status='Interesting Links';">Links</a> &nbsp;<font style="color: rgb(0, 80, 80);">|</font>&nbsp;""") 

            self.outmethod("""</td>""")
            self.outmethod("""<td><a href='http://grib.imim.es'><img src='%s/images/grib_logo.jpg' border=0></a></td>""" %(url))
            self.outmethod("""<td><a href='http://sbi.imim.es'><img src='%s/images/logo_sbi.gif' border=0></a></td>""" %(url))
            self.outmethod("""    </tr>""")
            self.outmethod("""    <tr><td class="links" align="center"></td></tr>""")
            self.outmethod("""  </tbody></table>""")

            self.outmethod("""  <table cellspacing=0 cellpadding=3 border=0 width=99%%>""")
            self.outmethod("""  <tr bgcolor="#eeedcd"><td></td></tr>""")
            self.outmethod("""  </table>""")

            # Print title of the report
            self.outmethod("<br><br><center><table border=0 style='border: 0;border-collapse: collapse;background-color:#CCCCFF;width: 600px;cellpadding=20'><tr><td style='padding:5px'><p style='text-align:center;font-weight: bold;color: #330066;font: 18px arial, sans-serif;'>%s</p></td></tr></table></center>\n" %(title))
            
            self.outmethod("<center><table border=0 style='border: 0;border-collapse: collapse;background-color:#CCCCCC;width: 600px;cellpadding=20'><tr><td style='padding:5px'>\n")


        return
        

    def closeReport(self, streamOutputMethod=sys.stdout.write, format="html"):
        """

        Closes the Biana report. This method is in charge of writing the closing HTML tags of the report

        "reportOutputMethod" is the write method where the report will be written (e.g. file_object.write)
              --> if no streamOutputMethod is given, report is written to stdout

        "reportFormat" is the formatting that will be applied to the report
             - 'html': report in html
             - 'txt': report in raw text

        Attention!!! This method does not close the file object associated to the report. The user is responsible
                     for closing the file (if needed).

        """
        if self.format == 'txt':
            self.outmethod("\n")

        elif self.format == "html":
            self.outmethod("</td></tr></table></center>\n  </body>\n</html>\n")

        # END OF elif self.format == "html":



        return 
        

    def addText(self, theText=None, type="regular" ):
        """
        This method adds text to the report.

        Depending on the 'type' of text, it will be written one way or another

          - 'regular' prints standard text
          - 'title' prints a section title

        """
        if self.format == 'txt':
            
            if type=="regular":
                self.outmethod("%s\n" %theText)
            elif type =="title":
                self.outmethod("%s\n\n" %(theText.upper()))

        elif self.format == "html":
            
            if type=="regular":
                self.outmethod("%s\n" %(theText))
            elif type =="title":
                self.outmethod("<h2>%s</h2>\n" %(theText))

        # END OF elif self.format == "html":



        return 
        

    def addResult(self, resultFileName=None, resultFilePath=None, associatedText=None, bulletedList=[] ):
        """
        This method adds information to the resport about a result file.

        'resultFileName' is the name you want to give to this result file

        'resultFilePath' is the path (or URL) of the result file (e.g. http://localhost/this_file.html

        'associatedText' is the text shown next to the results file name (i.e description of the file)

        'bulletedList' is a list of strings that will be shown under the results name as a bulleted list

        """
        if self.format == 'txt':
            self.outmethod("==> %s: %s\n" %(resultFileName, associatedText))
            self.outmethod("       - located in %s\n" %(resultFilePath))

            for one_associatedText in associatedText:
                self.outmethod("   - %s" %(associatedText))

        elif self.format == "html":
            self.outmethod("<li><p style='text-align:left;font-weight: bold;color: #000000;font: 14px arial, sans-serif;'><a href='%s'>%s</a>: %s</p></li>\n\n" %(resultFilePath, resultFileName,  associatedText))

            if bulletedList:
                self.outmethod("<ul>\n")
                for one_item in bulletedList:
                    self.outmethod("<li>%s</li>\n" %(one_item))
                self.outmethod("</ul>\n")

        # END OF elif self.format == "html":

        return 
    
