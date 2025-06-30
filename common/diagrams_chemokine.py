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

    def __init__(self, residue_list, protein_name, signal_start=None, signal_end=None, nobuttons=None):
    
        self.nobuttons = nobuttons if nobuttons is not None else 'arrestin'
        self.type = 'snakeplot'
        plot_data = {}
        plot_data['direction'] = [0, 0, 1, 0, 1, 0, 1, 0]  # 0: EC->IC, 1: IC->EC
        plot_data['helixRadius'] = 70

        self.receptorId = protein_name
        #self.family = protein_class
        self.output = ''

        # Store the signal range if provided.
        if signal_start is not None and signal_end is not None:
            print(signal_start, signal_end)
            self.signal_range = (signal_start, signal_end)
        else:
            self.signal_range = None


        self.sequence = residue_list
        self.segments = {}
        self.segments_full = OrderedDict()
        i = 0
        for r in self.sequence:
            if r.segment:
                segment = str(r.segment)
            elif r.segment_slug:  #from family aligment
                segment = str(r.segment_slug)

            if segment not in self.segments:
                self.segments[segment] = []
            self.segments_full[segment] = r.segment
            label = ''
            if r.ccn_number:
                label = str(r.ccn_number)
            elif hasattr(r, 'family_generic_number'):
                label = str(r.family_generic_number)

            displaylabel = f"{r.amino_acid}{r.sequence_number}"
            if r.ccn_number:
                displaylabel += f" {r.ccn_number}"
            if hasattr(r, 'frequency'):
                displaylabel += f"\n{r.frequency}"

            self.segments[segment].append([r.sequence_number, r.amino_acid, label, displaylabel])

            i += 1

        # Ensure Helix is present to avoid KeyError
        if "Helix" not in self.segments:
            self.segments["Helix"] = []

        # If C-term exists, append its residues to the Helix segment and remove it from segments.
        if "C-term" in self.segments:
            if "Helix" in self.segments:
                # Append C-term residues to the helix residues.
                self.segments["Helix"] += self.segments["C-term"]
            else:
                # If Helix doesn't exist, treat C-term as Helix.
                self.segments["Helix"] = self.segments["C-term"]
            # Remove the C-term so it doesn't get drawn as a terminal.
            del self.segments["C-term"]
            if "C-term" in self.segments_full:
                del self.segments_full["C-term"]


        for helix_num in range(1,2): #FIX for missing generic numbers
            rs = self.segments['Helix']
            for i in range(0,len(rs)):
                if not rs[i][2]:
                    if i+1<len(rs): #if there is a next one
                        if rs[i+1][2]: #if it has generic number
                            number = str(int(rs[i+1][2].split('x')[1])-1)
                            rs[i][2] = str(helix_num) + "x" + number

        self.helixWidth = 85            # Width of helix
        self.resNumPerRow = 4           # Residue number per row in helix
        self.angleDeg = 22.0            # Angle size of each helix turn
        self.residue_radius = 12        # Radius of the residue circle

        # svg image padding offset
        self.offsetX = 250 #-200
        self.offsetY = 0 #-50

        # margin between two helixes
        self.margin = 10

        # highest and lowest bound of this svg
        self.high = 0
        self.low = 0

        # keep track of max Y positions of intra/extra loops
        self.maxY = {'bottom':0,'top':0}
        self.maxX = {'left':0,'right':0}

        # helices length
        # helicesLength = Svg::getSnakePlotHelicesLength($baldwin, $helixWidth, $angleDeg) #FIXME

        # top and bottom residue coords in each helix
        self.TBCoords = {}

        self.output = ""
        self.traceoutput = ""
        self.helixoutput = ""

        # Draw sheets and helices
        self.count = 1
        self.count_sheet = 0
        for s in CHEMOKINE_SEGMENTS['Full']:
            if s in self.segments_full:

                if self.segments_full[s]=='Helix':
                    self.helixoutput += self.drawSnakePlotHelix(s)
                    self.count += 1
                if self.segments_full[s]=='B1':
                    self.helixoutput += self.drawSnakePlotSheet(s)
                    self.count += 1
                    self.count_sheet += 1
                if self.segments_full[s]=='B2':
                    self.helixoutput += self.drawSnakePlotSheet(s)
                    self.count += 1
                    self.count_sheet += 1
                if self.segments_full[s]=='B3':
                    self.helixoutput += self.drawSnakePlotSheet(s)
                    self.count += 1
                    self.count_sheet += 1
        
        # Draw loops
        self.count = 0
        for s in ['B1', '30s-loop', 'B2', '40s-loop', 'B3', '50s-loop', 'Helix']:
            if s in self.segments_full:
                if self.segments_full[s] == '30s-loop':
                    self.drawSnakePlotLoop(s)         
                elif self.segments_full[s] == '40s-loop':
                    self.drawSnakePlotLoop(s)
                elif self.segments_full[s] == '50s-loop':
                    self.drawSnakePlotLoop(s)
                else:
                    self.count += 1
            else:
                print(f"Warning: segment '{s}' not found in segments_full")

        self.drawCTermBox()
        
        # Draw terminals
        self.drawSnakePlotTerminals()



    def __str__(self):
        self.output = "<g id=snake transform='translate(0, " + str(-self.low+ self.offsetY) + ")'>" + self.traceoutput+self.output+self.helixoutput+self.drawToolTip() + "</g>"; #for resizing height
        return mark_safe(self.create(self.output,self.maxX['right']+30,self.high-self.low+self.offsetY*2,"snakeplot", self.nobuttons))


    def drawCTermBox(self):
        if not self.TBCoords:
            print("Warning: No helices drawn, skipping drawCTermBox.")
            return
        
        helix_num = max(self.TBCoords.keys())  # get the last helix drawn

        # Use the bottom coordinate of that helix as the last residue position.
        x, y = self.TBCoords[helix_num]['bottom']

        # Define an orientation/offset; here we use a positive offset (to the right)
        orientation = 1
        offset = 0

        # Calculate rectangle (box) position
        box_x = x + offset * orientation - 25  # center the box horizontally relative to the offset
        box_y = y + 30  # adjust vertically so it sits near the residue

        # Create the rectangle and text elements.
        box = (
            f"<rect class='C-term-box' x='{box_x}' y='{box_y}' rx='5' ry='5' "
            f"width='50' height='20' stroke='black' fill='white' stroke-width='1'/>"
        )
        text = (
            f"<text class='C-term-box' x='{x + offset * orientation}' y='{y+43}' text-anchor='middle' "
            f"font-size='12' font-family='helvetica'>C-term</text>"
        )

        # Append them to the SVG output.
        self.output += box + text

    def drawSnakePlotTerminals(self):
        y_offset = 50
        font_size = 12
        font_family = 'helvetica'
        bezier_pull = 80
        between_residues = 18

        for i in ['N', 'C']:
            drawn_residues = []
            name = i + "-term"
            if name not in self.segments:
                continue  # continue if no terminus

            rs = self.segments[name]  # get residues

            # TEMP FIX for N-term segments to concatenate
            if name == "N-term":
                rs = self.segments["N-term"] + self.segments["CX"] + self.segments["N-loop"]

            if i == 'N':
                orientation = 1
                y_max = self.maxY['bottom'] - between_residues * 4
                x_max = self.maxX['left']
                position = 'bottom'
                linked_helix = 1
                rs.reverse()
            else:
                orientation = 1
                y_max = self.maxY['bottom'] + between_residues * 4
                x_max = self.maxX['left']
                position = 'bottom'
                linked_helix = 4

            x1 = self.TBCoords[linked_helix][position][0]
            y1 = self.TBCoords[linked_helix][position][1]

            # For N-term, use the original parameters:
            if name == "N-term":
                x2 = x1 - 90 * orientation
                y2 = y_max + 90
                bezierX = x1 + 60 * orientation
                bezierY = (y_max + y1) / 2 + 60 * orientation
            
            else:  # C-term
                x2 = x1 - 90 * orientation
                y2 = y_max + 90
                bezierX = x1 + 60 * orientation
                bezierY = (y_max + y1) / 2 + 60 * orientation

                # Filter residues to only those whose segment is "C-term"
                c_term_residues = [r for r in rs if str(getattr(r, 'segment', "")) == "C-term"]
                # Only adjust (scale) the curve if there are less than 12 C-term residues.
                if len(c_term_residues) < 1:
                    current_length = self.lengthbezier([x1, y1], [bezierX, bezierY], [x2, y2], 0.001)
                    target_length = len(rs) * between_residues + 30
                    scale_factor = target_length / float(current_length)
                    x2 = x1 + (x2 - x1) * scale_factor
                    y2 = y1 + (y2 - y1) * scale_factor
                    bezierX = x1 + (bezierX - x1) * scale_factor
                    bezierY = y1 + (bezierY - y1) * scale_factor

            points = "M " + str(x1) + " " + str(y1) + " Q" + str(bezierX) + " " + str(bezierY) + " " + str(x2) + " " + str(y2)
            # Draw the bezier path for terminal connection:
            self.output += "<path class='" + name + " long' d='" + points + "' stroke='black' fill='none' stroke-width='2' />"
            pos = 40
            length = self.lengthbezier([x1, y1], [bezierX, bezierY], [x2, y2], 0.001)

            bend = 0
            distance_between_rows = 30
            pos_bend = 0
            bend_direction = -1 * orientation

            for j in range(0, len(rs)):
                r = rs[j]
                if pos < length:
                    where = self.wherebezier([x1, y1], [bezierX, bezierY], [x2, y2], 0.001, pos)
                else:
                    if pos_bend == 0 and bend != 0:
                        where[1][0] = where[1][0] - between_residues * bend_direction
                        where[1][1] = where[1][1] + orientation * distance_between_rows / 2
                    elif pos_bend == between_residues and bend != 0:
                        where[1][0] = where[1][0] + between_residues * bend_direction
                        where[1][1] = where[1][1] + orientation * distance_between_rows / 2
                    else:
                        where[1][0] = where[1][0] + between_residues * bend_direction
                        where[1][1] = where[1][1]
                    pos_bend += between_residues
                    if pos_bend >= abs(x2 - x_max) - 40:
                        pos_bend = 0
                        bend += 1
                        bend_direction = -bend_direction

                if bend == 0:
                    labely = where[1][1]

                display_label = rs[j][2] if name == "N-term" else rs[j][3]
                drawn_residues.append(
                    self.DrawResidue(
                        x=where[1][0],                 # X coordinate
                        y=where[1][1],                 # Y coordinate
                        aa=r[1],                       # Amino acid letter
                        residue_number=r[0],           # Sequence number
                        label=display_label,           # Tooltip label (should be the CCN number or other info)
                        radius=self.residue_radius - 1,# Circle radius
                        resclass=name + " long"        # SVG class
                    )
                )
                pos += between_residues

                if where[1][1] < self.low:
                    self.low = where[1][1]
                if where[1][1] > self.high:
                    self.high = where[1][1]

            if name == 'N-term':
                drawn_residues = drawn_residues[::-1]
            self.output += ''.join(drawn_residues)
            # Draw terminal label without any onclick toggling:
            self.output += "<rect class='" + name + " long segment' x='" + str(self.TBCoords[linked_helix][position][0] + 50 * orientation - 25) + "' y='" + str((labely + self.TBCoords[linked_helix][position][1]) / 2 - 13) + "' rx='5' ry='5' width='50' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
            self.output += "<text class='" + name + " long segment' x='" + str(self.TBCoords[linked_helix][position][0] + 50 * orientation) + "' y='" + str((labely + self.TBCoords[linked_helix][position][1]) / 2) + "' text-anchor='middle' font-size='" + str(font_size) + "' font-family='" + font_family + "'>" + name + "</text>"


    def drawSnakePlotHelix(self, segment):
        rs = self.segments[segment]
        helix_num = self.count
        self.TBCoords[helix_num] = {}

        if helix_num%2!=0: rs.reverse() # reverse direction for even helix because they go from inside to outside

        output_residues = []

        res_num = len(self.segments[segment])
        output_residue_in = ''
        output_residue_out = ''
        output_trace = ''

        startX = self.helixWidth+self.offsetX+(self.margin+self.helixWidth)*(helix_num-1)-(self.count_sheet*20)
        startY = self.offsetY

        row_length = 3
        row_pos = 0
        row = 0
        prevGeneric = 0
        bulgeX = 0
        bulgeY = 0
        bulge = 0
        skip = 0
        indentX = -self.residue_radius+3
        indentY = 3
        for i in range(0,res_num):

            # move left as you go down a row
            x = round(startX-row_pos*self.residue_radius*1.6+indentX+bulgeX)

            # Move down with right amount
            y = round(startY+row*self.residue_radius*2.4+row_pos*self.residue_radius*0.5+indentY+bulgeY)
            output_residue = self.DrawResidue(x,y,rs[i][1], rs[i][0], rs[i][2], self.residue_radius)


            if x<self.maxX['left']: self.maxX['left'] = x
            if x>self.maxX['right']: self.maxX['right'] = x

            row_pos += 1
            if bulge==1:
                if row_pos==1:  # if first in row, use space for bulge
                    bulgeY = -3
                    bulgeX = 10
                else:
                    bulgeY = -3
                    bulgeX = 7
                rs[i][2] = prevGeneric # make it the prev one, to catch missing ones correctly
                bulge = 0

            if row_length==3:
                output_residue_in += output_residue
            else:
                output_residue_out += output_residue

            output_residues.append(output_residue)

            if i==0: self.TBCoords[helix_num]['top'] = [x,y]
            if i==res_num-1: self.TBCoords[helix_num]['bottom'] = [x,y]


            if (row_pos==1 and row!=0) or (skip==1 and row_pos==2): # if need for trace
                if row_length==3: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-21)+" "+str(y-8)+" T"+str(x)+" "+str(y)
                if row_length>=4: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-24)+" "+str(y-7)+" T"+str(x)+" "+str(y)
                output_trace += "<path d='" + points + "' stroke='grey' fill='none' stroke-width='2'  />"

            # alternate between 4 and 3 res per row
            if row_length>3 and row_pos>=row_length:
                row_length=3
                row_pos = 0
                row += 1
                bulgeX = 0
                bulgeY = 0
                indentX = -self.residue_radius+3
                indentY = 3
            elif row_length==3 and row_pos>=3:
                row_length=4
                row_pos = 0
                row += 1
                bulgeX = 0
                bulgeY = 0
                indentX = 0
                indentY = 0

            skip = 0
            prevX = x
            prevY = y
            prevGeneric = rs[i][2]

        temp = ''
        if helix_num%2!=0: output_residues.reverse()
        for res in output_residues:
            temp += res

        return output_trace+temp

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

            output_residue = self.DrawResidueSquare(x, y, rs[i][1], rs[i][0], rs[i][2], self.residue_radius)
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

            self.output += self.DrawResidue(where[1][0],where[1][1],r[1], r[0], r[2], self.residue_radius-1,name)
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