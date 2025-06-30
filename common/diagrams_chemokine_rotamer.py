"""
Modified Script for Chemokine SVG Diagram

Original Source:
    Title: diagrams_gprotein.py
    Author: University of Copenhagen
    Date Accessed: September 2023
    Original License: Apache License 2.0
    Original Copyright: © 2015 University of Copenhagen
    Original Source Location: https://github.com/protwis/protwis/blob/master/common/diagrams_gprotein.py

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
    either express or implied. See the License for the specific language governing permissions and limitations under the License.

Description of Modifications:
    - Adapted the original G protein plot style to represent chemokines in `DrawArrestinPlot`.
    - Updated segment data with `CHEMOKINE_SEGMENTS` for chemokine-specific regions.
    - Modified `drawSnakePlotHelix`, `drawSnakePlotSheet`, `drawSnakePlotLoop`, and `drawSnakePlotTerminals` to reflect chemokine structure.
    - Adjusted SVG output functions to better fit chemokine visual representation.
    - Creating a separate `DrawArrestinPlot` class focused on chemokine data.

NOTICE: This file includes derivative works of the original code licensed under the Apache License 2.0 by the University of Copenhagen. 
Modifications made to the original code are marked by the description above. This modified version retains the same license.
"""


from common.diagrams import Diagram
#from common.definitions import CHEMOKINE_SEGMENTS

from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from django.utils.safestring import mark_safe

from math import cos, sin, pi, floor, sqrt
from datetime import datetime
from collections import OrderedDict

CHEMOKINE_SEGMENTS = OrderedDict([
        ('Full', ['N-term', 'CX', 'N-loop', 'B1', '30s-loop', 'B2',
       '40s-loop', 'B3', '50s-loop', 'Helix', 'C-term']),
        ('Structured', ['B1', 'B2', 'B3',
       'Helix']),
    ])

# GPCRdb snakeplot style modified for chemokines
# https://github.com/protwis/protwis
class DrawArrestinPlot(Diagram):
    def __init__(self, rotamer_list, protein_name, nobuttons=None):
        self.nobuttons = 'arrestin'
        self.type = 'snakeplot'
        self.receptorId = protein_name
        self.output = ''
        self.sequence = rotamer_list  # Use rotamers instead of residues
        self.segments = {}
        self.segments_full = OrderedDict()

        # Parse rotamers into segments
        for r in self.sequence:
            segment = r.segment or "Unknown"
            if segment not in self.segments:
                self.segments[segment] = []
            self.segments_full[segment] = r.segment
            label = r.generic_number or ""
            displaylabel = f"{r.amino_acid}{r.sequence_number} \n {label}"
            self.segments[segment].append([r.sequence_number, r.amino_acid, label, displaylabel])

        # Ensure generic numbers for helices
        if 'Helix' in self.segments:
            rs = self.segments['Helix']
            for i in range(len(rs)):
                if not rs[i][2] and i + 1 < len(rs) and rs[i + 1][2]:
                    number = int(rs[i + 1][2].split('x')[1]) - 1
                    rs[i][2] = f"1x{number}"

        # SVG drawing configuration
        self.helixWidth = 85
        self.resNumPerRow = 4
        self.angleDeg = 22.0
        self.residue_radius = 12
        self.offsetX = 250
        self.offsetY = 0
        self.margin = 10
        self.high = 0
        self.low = 0
        self.maxY = {'bottom': 0, 'top': 0}
        self.maxX = {'left': 0, 'right': 0}
        self.TBCoords = {}
        self.traceoutput = ""
        self.helixoutput = ""

        # Draw SVG components
        self.count = 1
        self.count_sheet = 0
        for segment in CHEMOKINE_SEGMENTS['Full']:
            if segment in self.segments_full:
                if self.segments_full[segment] == 'Helix':
                    self.helixoutput += self.drawSnakePlotHelix(segment)
                    self.count += 1
                elif self.segments_full[segment] in ['B1', 'B2', 'B3']:
                    self.helixoutput += self.drawSnakePlotSheet(segment)
                    self.count += 1
                    self.count_sheet += 1

        # Draw loops and terminals
        self.count = 0
        for segment in ['B1', '30s-loop', '40s-loop', '50s-loop', 'Helix']:
            if segment in self.segments_full:
                self.drawSnakePlotLoop(segment)
        self.drawSnakePlotTerminals()

    def __str__(self):
        self.output = (
            f"<g id=snake transform='translate(0, {-self.low + self.offsetY})'>"
            f"{self.traceoutput}{self.output}{self.helixoutput}{self.drawToolTip()}</g>"
        )
        return mark_safe(self.create(self.output, self.maxX['right'] + 30, self.high - self.low + self.offsetY * 2, "snakeplot", self.nobuttons))



    def drawSnakePlotTerminals(self):
        y_offset = 50
        font_size = 12
        font_family = 'helvetica'
        bezier_pull = 80

        between_residues = 18


        for i in ['N','C']:
            drawn_residues = []

            name = i+"-term"
            if name not in self.segments: continue # continue if no terminus

            rs = self.segments[name] # get residues

            ### TEMP FIX for N-term segments to concatenate
            if name=="N-term":
                rs = self.segments["N-term"] + self.segments["CX"] + self.segments["N-loop"]
            ###

            if i=='N':
                orientation = 1
                y_max = self.maxY['bottom']-between_residues*4
                x_max = self.maxX['left']
                position = 'bottom'
                linked_helix = 1
                rs.reverse()
            else:
                orientation = 1
                y_max = self.maxY['bottom']+between_residues*4
                x_max = self.maxX['left']
                position = 'bottom'
                linked_helix = 4

            x1 = self.TBCoords[linked_helix][position][0]
            y1 = self.TBCoords[linked_helix][position][1]

            # Get positions of two  linking residues from each helix
            x2 = x1-30
            y2 = y1+60*orientation

            # Make line and box for short version
            points = "M "+str(x1)+" "+str(y1)+" Q"+str(x1+30)+" "+str(y2)+" "+str(x2)+" "+str(y2)
            self.output += "<path class='"+name+" short' d='" + points + "' stroke='black' fill='none' stroke-width='2' />"
            self.output += "<rect class='"+name+" short segment' onclick='toggleLoop(\"."+name+"\",\"short\",false,this);' x="+str(x2-25)+" y="+str(y2-13)+" rx=5 ry=5 width='50' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
            self.output += str("<text class='"+name+" short segment' onclick='toggleLoop(\"."+name+"\",\"short\",false,this);' x="+str(x2)+" y="+str(y2)+" text-anchor='middle' font-size="+str(font_size)+" font-family='"+font_family+"'>"+name+"</text>")

            x2 = x1-90*orientation
            y2 = y_max+90
            bezierX = x1+60*orientation
            bezierY = (y_max+y1)/2+60*orientation

            points = "M "+str(x1)+" "+str(y1)+" Q"+str(bezierX)+" "+str(bezierY)+" "+str(x2)+" "+str(y2)

            pos = 40

            length = self.lengthbezier([x1,y1],[bezierX,bezierY],[x2,y2],0.001)

            bend = 0
            distance_between_rows = 30
            pos_bend = 0
            bend_direction = -1*orientation

            for i in range(0,len(rs)):

                r = rs[i]
                if pos<length:
                    where = self.wherebezier([x1,y1],[bezierX,bezierY],[x2,y2],0.001,pos)
                else:
                    if pos_bend==0 and bend!=0: #if first residue in line put in middle
                        where[1][0] = where[1][0]-between_residues*bend_direction
                        #where[1][0] = where[1][0]
                        where[1][1] = where[1][1]+orientation*distance_between_rows/2
                    elif pos_bend==between_residues and bend!=0: #if 2nd residue in line put in middle
                         #where[1][0] = where[1][0]-between_residues*bend_direction
                         where[1][0] = where[1][0]+between_residues*bend_direction
                         where[1][1] = where[1][1]+orientation*distance_between_rows/2
                    else:
                        where[1][0] = where[1][0]+between_residues*bend_direction
                        where[1][1] =  where[1][1]
                    last_bend_x = where[1][0]
                    last_bend_y = where[1][1]

                    pos_bend += between_residues
                    if pos_bend>=abs(x2-x_max)-40: #no more bend left
                        pos_bend = 0
                        bend += 1
                        if bend_direction==1:
                            bend_direction = -1
                        elif bend_direction==-1:
                            bend_direction = 1

                if i==0: self.output += "<line class='"+name+" long' x1="+str(x1)+" y1="+str(y1)+" x2="+str(where[1][0])+" y2="+str(where[1][1])+" stroke='black' fill='none' stroke-width='2' stroke-dasharray2='1,1' />"

                if bend==0: labely = where[1][1]

                drawn_residues.append(self.DrawResidue(where[1][0],where[1][1],r[1], r[0], rs[i][3], self.residue_radius-1,name+" long"))
                pos += between_residues

                if where[1][1]<self.low: self.low = where[1][1]
                if where[1][1]>self.high: self.high = where[1][1]

            if name=='N-term': drawn_residues = drawn_residues[::-1]
            self.output += ''.join(drawn_residues)
            self.output += "<rect onclick='toggleLoop(\"."+name+"\",\"long\",false,this);' class='"+name+" long segment' x="+str(self.TBCoords[linked_helix][position][0]+50*orientation-25)+" y="+str((labely+self.TBCoords[linked_helix][position][1])/2-13)+" rx=5 ry=5 width='50' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
            self.output += str("<text onclick='toggleLoop(\"."+name+"\",\"long\",false,this);' class='"+name+" long segment' x="+str(self.TBCoords[linked_helix][position][0]+50*orientation)+" y="+str((labely+self.TBCoords[linked_helix][position][1])/2)+" text-anchor='middle' font-size="+str(font_size)+" font-family='"+font_family+"'>"+name+"</text>")



    def drawSnakePlotHelix(self, segment):
        rs = self.segments[segment]
        helix_num = self.count
        self.TBCoords[helix_num] = {}
        if helix_num % 2 != 0:
            rs.reverse()

        output_residues = []
        startX = self.helixWidth + self.offsetX + (self.margin + self.helixWidth) * (helix_num - 1) - (self.count_sheet * 20)
        startY = self.offsetY
        row_length = 3
        row_pos = 0
        row = 0

        for i, residue in enumerate(rs):
            x = round(startX - row_pos * self.residue_radius * 1.6 - self.residue_radius + 3)
            y = round(startY + row * self.residue_radius * 2.4 + row_pos * self.residue_radius * 0.5 + 3)
            output_residue = self.DrawResidue(x, y, residue[1], residue[0], residue[3], self.residue_radius)
            output_residues.append(output_residue)

            if i == 0:
                self.TBCoords[helix_num]['top'] = [x, y]
            if i == len(rs) - 1:
                self.TBCoords[helix_num]['bottom'] = [x, y]

            row_pos += 1
            if row_length > 3 and row_pos >= row_length:
                row_length = 3
                row_pos = 0
                row += 1
            elif row_length == 3 and row_pos >= 3:
                row_length = 4
                row_pos = 0
                row += 1

        return ''.join(output_residues)

    def drawSnakePlotSheet(self, segment):
        rs = self.segments[segment]
        helix_num = self.count
        self.TBCoords[helix_num] = {}

        if helix_num % 2 != 0: 
            rs.reverse()  # Reverse direction for even helix

        output_residues = []
        res_num = len(self.segments[segment])
        startX = 50 + self.offsetX + (self.margin + self.helixWidth) * (helix_num - 1) - (self.count_sheet * 20)
        startY = self.offsetY

        for i in range(0, res_num):
            # Calculate the position of each residue in the sheet
            x = round(startX)
            y = round(startY + i * self.residue_radius * 1.5)

            output_residue = self.DrawResidueSquare(x, y, rs[i][1], rs[i][0], rs[i][3], self.residue_radius)
            output_residues.append(output_residue)

            if x < self.maxX['left']: 
                self.maxX['left'] = x
            if x > self.maxX['right']: 
                self.maxX['right'] = x

            if i == 0: 
                self.TBCoords[helix_num]['top'] = [x, y]
            if i == res_num - 1: 
                self.TBCoords[helix_num]['bottom'] = [x, y]

        return ''.join(output_residues)


    def drawSnakePlotLoop(self, segment):
        y_offset = 50
        font_size = 12
        font_family = 'courier'
        bezier_pull = 80
        name = segment
        x_at_max_y = 0

        rs = self.segments[segment] # get residues

        start = 1
        res_before = []
        res_helix = []
        res_after = []

        if self.count % 2 == 0:
            position = 'bottom'
            orientation = 1
        if segment == "40s-loop":
            position = 'bottom'
            orientation = 1
        else:
            position = 'top'
            orientation = -1

        if self.count not in self.TBCoords:
            return 0

        if self.count+1 not in self.TBCoords:
            return 0

        # Get positions of two  linking residues from each helix
        x1 = self.TBCoords[self.count][position][0]
        y1 = self.TBCoords[self.count][position][1]
        x2 = self.TBCoords[self.count+1][position][0]
        y2 = self.TBCoords[self.count+1][position][1]

        boxX = (x1+x2)/2 # midway between
        if position=='top':
            boxY = min(y1,y2)-y_offset # over helix
            y_indent = -1*bezier_pull
        if position=='bottom':
            boxY = max(y1,y2)+y_offset # over helix
            y_indent = bezier_pull

        points = str(x1)+","+str(y1)+" "+str(boxX)+","+str(boxY)+" "+str(x2)+","+str(y2)
        points2 = "M "+str(x1)+" "+str(y1)+" Q"+str(boxX)+" "+str(boxY+y_indent)+" "+str(x2)+" "+str(y2)

        # Getting midpoint of Bezier curve http://www.svgbasics.com/curves.html
        Dx = ((x1+boxX)/2)
        Ex = ((x2+boxX)/2)
        Fx = (Dx+Ex)/2

        Dy = ((y1+boxY+y_indent)/2)
        Ey = ((y2+boxY+y_indent)/2)
        Fy = (Dy+Ey)/2

        y_indent = y_indent*len(rs)/5 # get an approx need for y_indent for size of loop

        loop_long_length = 0
        super_loop_long_length = 40
        between_residues = 18

        length_of_residues_in_loop = len(rs)*between_residues-self.residue_radius
        length = self.lengthbezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001)

        if len(rs)<super_loop_long_length:
            tries = 0 # adjust size
            while abs(length-length_of_residues_in_loop-70)>5:
                # print(abs(length-length_of_residues_in_loop+100),length,length_of_residues_in_loop,tries)
                if length-length_of_residues_in_loop-70>5:
                    y_indent *=0.9
                else:
                    y_indent *=1.1
                length = self.lengthbezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001)

                tries += 1
                if tries>100:
                    break

        pos = (length-length_of_residues_in_loop)/2 # get start pos

        indentX = 0
        indentY2 = 0
        prev_where = [x1,y1]

        # make rounded arc
        points2 = "M "+str(x1)+" "+str(y1)+" Q"+str(boxX)+" "+str(boxY+y_indent)+" "+str(x2)+" "+str(y2)
        labelbox = self.wherebezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001,length/2)

        labelbox[1][1] += orientation*40

        self.output += "<path class='"+name+"' d='" + points2 + "' stroke='black' fill='none' stroke-width='2' />"

        max_y = y1
        for i in range(0,len(rs)):
            r = rs[i]
            where = self.wherebezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001,pos)

            self.output += self.DrawResidue(where[1][0],where[1][1],r[1], r[0], r[3], self.residue_radius-1,name)
            pos += between_residues

            if where[1][1]>self.high: self.high = where[1][1]
            if where[1][1]<self.low: self.low = where[1][1]
            prev_where = where[1][0],where[1][1]

            if orientation==-1:
                if where[1][1]<self.maxY[position]: self.maxY[position] = where[1][1]
            else:
                if where[1][1]>self.maxY[position]: self.maxY[position] = where[1][1]

            if orientation==-1:
                if where[1][1]<max_y:
                    max_y = where[1][1]
                    x_at_max_y = where[1][0]
            else:
                if where[1][1]>max_y:
                    max_y = where[1][1]
                    x_at_max_y = where[1][0]
            x_at_max_y = where[1][0]

        if orientation==1:
            max_y = max_y+25
        else:
            max_y = max_y-20
        self.output += "<rect onclick='toggleLoop(\"."+name+"\",\"long\");' class='"+name+"' x="+str(x_at_max_y-33)+" y="+str(max_y-13)+" rx=5 ry=5 width='65' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
        self.output += str("<text  onclick='toggleLoop(\"."+name+"\",\"long\");' class='"+name+"' x="+str(x_at_max_y)+" y="+str(max_y)+" text-anchor='middle' font-size="+str(font_size)+" font-family='"+font_family+"'>"+name+"</text>")
