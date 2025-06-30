"""
Modified Script for SVG Diagram Generation with Chemokine Plotting

Original Source:
    Title: diagrams_gpcr.py
    Author: University of Copenhagen
    Date Accessed: September 2023
    Original License: Apache License 2.0
    Original Copyright: © 2015 University of Copenhagen
    Original Source Location: https://github.com/protwis/protwis/blob/master/common/diagrams.py

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
    either express or implied. See the License for the specific language governing permissions and limitations under the License.

Description of Modifications:
    - Refactored core SVG creation methods adapted for chemokine plotting.
    - Consolidated and enhanced `drawColorPanel`, `create`, `DrawResidue`, and `DrawBackbone`.
    - Modified `DrawGproteinPlot` to separate G protein-specific functions like `drawSnakePlotHelix`, `drawSnakePlotLoop`, and `drawSnakePlotSheet`.

NOTICE: This file includes derivative works of the original code licensed under the Apache License 2.0 by the University of Copenhagen. 
Modifications made to the original code are marked by the description above. This modified version retains the same license.
"""


from math import cos, sin, tan, pi, sqrt, pow
import string, time, math, random


def uniqid(prefix='', more_entropy=False):
    m = time.time()
    uniqid = '%8x%05x' %(int(math.floor(m)),int((m-math.floor(m))*1000000))
    if more_entropy:
        valid_chars = list(set(string.hexdigits.lower()))
        entropy_string = ''
        for i in range(0,10,1):
            entropy_string += random.choice(valid_chars)
        uniqid = uniqid + entropy_string
    uniqid = prefix + uniqid
    return uniqid

class Diagram:
    def create(self, content, sizex, sizey, name, nobuttons):
        #diagram_js = self.diagramJS()
        if nobuttons == 'chemokine':
            return (
                "<svg id=\"" + str(name) + "\" " +
                "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" + str(sizex) + "\" height=\"" + str(sizey) + "\" " +
                "style='stroke-width: 0px; background-color: white;'>\n" + content + "</svg>" +
                self.drawColorPanel()
            )
        elif nobuttons:
            return (
                "<svg id=\"" + str(name) + "\" " +
                "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" + str(sizex) + "\" height=\"" + str(sizey) + "\" " +
                "style='stroke-width: 0px; background-color: white;'>\n" + content + "</svg>" +
                self.drawColorPanel(nobuttons)
            )
        else:
            return (
                "<svg id=\"" + str(name) + "\" " +
                "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" + str(sizex) + "\" height=\"" + str(sizey) + "\" " +
                "style='stroke-width: 0px; background-color: white;'>\n" + content + "</svg>" +
                self.drawColorPanel()
            )

    def drawToolTip(self):
        output2 = """<g id='tool-tip-{}' transform='translate(0,0)' visibility='hidden'>
            <rect x='0' y='-40' width='1' height='25' stroke='black' fill='white' stroke-width='1' />
            <text x='0' y='-23' text-anchor='middle' font-family='Arial' font-size='12' fill='black'></text>
            </g>""".format(self.type)
        output= ""

        return output

    def drawColorPanel(self, nobuttons=None):
    
        boxstyle = """<style>
        .pick-color  {
          display:inline-block;
          width: 35px;
          height: 20px;
          margin: 1px;
          border-radius: 5px;
          border: 2px solid #000;
        }
    
        .tooltip-inner {
            white-space:pre-wrap;
            max-width:none;
        }
        </style>
        """
        
        # If nobuttons is True or truthy, just return style (or even an empty string if you want nothing)
        if nobuttons:
            return boxstyle
    
        presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                                    'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                                    'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                                    'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                                    'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                                    'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                                    'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                                    'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                                    'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                                    'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                                    '-': ['#FFFFFF', '#000000']
                                    }
        fillcolors = [['#CCCCCC', '#000000']]
        for key, value in presetColors.items():
            if value not in fillcolors:
                fillcolors.append(value)
    
        colors = "<div id=\"cp2_" + self.type + "\" class=\"\">"
        for color in fillcolors:
            colors += "<div class='pick-color " + self.type + "' id='pick-" + color[0] + "-" + color[1] + "' style='background-color: " + color[0] + ";'>&nbsp;</div>"
    
        colors += "<div class=\"input-group-addon\" style='width:0;padding:0;border:0;background-color:0;display:inline;position:relative;top:-6px'><i class='pick-color " + self.type + " selected'></i></div>"
        colors += "<input type=\"text\" id='custom_color_" + self.type + "' value=\"#00AABB\" class=\"\" size=8 />"
        colors += "</div>"
    
        output = ("<br>Pick color:" + colors)
    
        # Append the Properties and Clear buttons plus our new Signal Sequence button.
        output += ('<br><button style="width:120px;" onclick="applyPresentColors(\'' + self.type + '\')">Properties</button> ' +
                   '<button style="width:120px;" onclick="resetColors(\'' + self.type + '\')">Clear</button> ' +
                   '<button style="width:120px;" onclick="toggleSignalSequence(\'' + self.type + '\')">Signal Sequence</button>')
    
        return boxstyle + output


    def DrawResidue(self, x, y, aa, residue_number, label, radius, resclass='', cfill="white", precolor=False, signal=False):
        """
        Draws a circular SVG residue with the generic residue number as the id.
        If `signal` is True, adds a data-signal attribute.

        For residues in the N-terminal (as detected by 'N-term' in resclass) and if a signal_range is defined,
        this method will automatically set the signal flag to True if the residue number is within the range.
        """
        generic_id = label.replace('.', '_').replace('x', '_')
        
        # Check if this residue is part of the N-terminal and if a signal_range is defined.
        if "N-term" in resclass and self.signal_range:
            try:
                # Ensure residue_number and the range boundaries are compared as integers.
                if  int(self.signal_range[0]) <= int(residue_number) <= int(self.signal_range[1]):
                    signal = True
            except ValueError:
                # If residue_number cannot be converted to an int, ignore the check.
                pass

        idtext = f"{generic_id}t"
        x = round(x)
        y = round(y)
        # If signal is True, include the custom attribute; otherwise leave it empty.
        data_signal = 'data-signal="true"' if signal else ''
        
        # ---- TOOLTIP UPDATE ----
        tooltip = f"Seq: {residue_number}\nCCN: {label}"
        
        output = f"""
            <circle class='{resclass} rcircle' cx='{x}' cy='{y}' r='{radius}' stroke='black' stroke-width='2' fill='{cfill}'
            fill-opacity='1' id='{generic_id}' title='{tooltip}' original_title='{label}' {data_signal} />
            <text x='{x}' y='{y+6}' text-anchor='middle' font-family='helvetica' font-size='16' fill=''
            id='{idtext}' class='rtext {resclass}' title='{tooltip}' original_title='{label}' {data_signal}> {aa} </text>
        """
        return output




    def DrawResidueSquare(self, x, y, aa, residue_number, label, radius, resclass='', cfill="white", precolor=False):
        """
        Draws a square SVG residue with the generic residue number as the id.
        """
        generic_id = label.replace('.', '_').replace('x', '_')  # Replace special characters
        idtext = f"{generic_id}t"

        x = round(x)
        y = round(y)
        
        # ---- TOOLTIP UPDATE ----
        tooltip = f"Seq: {residue_number}\nCCN: {label}"

        output = f"""
            <rect class='{resclass} rcircle' x='{x-radius*(1.5/2)}' y='{y-radius*(1.5/2)}' height='{radius*1.5}' width='{radius*1.5}' 
            stroke='black' stroke-width='1' fill='{cfill}' fill-opacity='1' id='{generic_id}' title='{tooltip}' original_title='{label}' />
            <text x='{x}' y='{y+6}' text-anchor='middle' font-family='helvetica' font-size='16' fill=''
            id='{idtext}' class='rtext {resclass}' title='{tooltip}' original_title='{label}'> {aa} </text>
        """
        return output

    
    def deg2rad(self,degrees):
        radians = pi * degrees / 180
        return radians

    def bezier(self,p0,p1,p2,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        v1x = p1[0]-p0[0]
        v1y = p1[1]-p0[1]

        i1 = [p0[0]+(p1[0]-p0[0])*t,p0[1]+(p1[1]-p0[1])*t]
        i2 = [p1[0]+(p2[0]-p1[0])*t,p1[1]+(p2[1]-p1[1])*t]

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def bezier_high(self,p0,p1,p2,p3,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        i1 = self.bezier(p0,p1,p2,t)
        i2 = self.bezier(p1,p2,p3,t)

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def bezier_high2(self,p0,p1,p2,p3,p4,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        #https://www.jasondavies.com/animated-bezier/
        i1 = self.bezier_high(p0,p1,p2,p3,t)
        i2 = self.bezier_high(p1,p2,p3,p4,t)

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def lengthbezier(self,p0,p1,p2,step,p3=False,p4=False):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        pos = 0
        length = 0
        p = p0
        while pos <= 1:


            if p3==False:
                xy = self.bezier(p0,p1,p2,pos)
            elif p4==False:
                xy = self.bezier_high(p0,p1,p2,p3,pos)
            elif p4!=False:
                xy = self.bezier_high2(p0,p1,p2,p3,p4,pos)

            length += math.sqrt( (xy[0]-p[0])**2 + (xy[1]-p[1])**2 )
            p = xy
            pos += step

        return round(length)

    def wherebezier(self,p0,p1,p2,step,stop,p3=False,p4=False):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        pos = 0
        length = 0
        p = p0
        xy = [0,0]

        if stop<0:
            if p3==False:
                stop = self.lengthbezier(p0,p1,p2,step)+stop
            elif p4==False:
                stop = self.lengthbezier(p0,p1,p2,step,p3)+stop
            else:
                stop = self.lengthbezier(p0,p1,p2,step,p3)+stop

        while pos <= 1:

            if length>stop: #stop if it reached the length along the line
                break

            if p3==False:
                xy = self.bezier(p0,p1,p2,pos)
            elif p4==False:
                xy = self.bezier_high(p0,p1,p2,p3,pos)
            else:
                xy = self.bezier_high2(p0,p1,p2,p3,p4,pos)

            length += math.sqrt( (xy[0]-p[0])**2 + (xy[1]-p[1])**2 )
            p = xy
            pos += step

        return pos,xy


    #find slope and y-intercept of a line through two points
    def LineEquation(self,p1, p2):
        # slope
        m = (p2['y']-p1['y'])/(p2['x']-p1['x'])

        # y-intercept
        b = p1['y'] - m*p1['x']

        # direction
        if p1['y'] >= p2['y'] and p1["x"] > p2["x"]:
            x = -1
            y = -1
        elif p1['y'] > p2['y'] and p1["x"] < p2["x"]:
            x = 1
            y = 1
        elif p1['y'] < p2['y'] and p1["x"] > p2["x"]:
            x = -1
            y = -1
        elif p1['y'] <= p2['y'] and p1["x"] < p2["x"]:
            x = 1
            y = 1

        return {'m':m, 'b':b, 'x':x, 'y':y}

    def MoveAlongLine(self,pixels_to_move,m,perpendicular, xDir=1, yDir=1):
        # perpendicular line
        if perpendicular == True:
            if m != 0:
                m = -1/m
            else:
                return {"x":0, "y":pixels_to_move}

        # a^2+b^2=c^2 where a=x and b=mx
        # x^2+mx^2=c^2 (the b of y=mx+b can be omitted as the y-intercept is zero)
        c = pixels_to_move
        x = sqrt( pow(c,2)/(1+pow(m,2)) )
        y = m*x

        return {"x":x*xDir, "y":y*yDir}

    #Draws the full backbone representation for one (box) helix

    def DrawBackbone(self,coordinates):

        # box properties
        numSides = 4

        # start at residue position 1
        cur_res = 1

        # loop through residues and draw backbone
        output = ""
        #for (i=1;i<=count(coordinates);i++) {
        for i in range(1,len(coordinates)+1):
            cur_res = i

            # find next residue
            next_res = cur_res-1
            if next_res<1:
                next_res = 19

            # find prev residue
            prev_res = cur_res+1;
            if prev_res>len(coordinates):
                prev_res = 2

            # next-in-wheel residue orientation
            wheel_next_res = cur_res-numSides
            if wheel_next_res <= (1-numSides):  # if current res is the last res, next is the first res
                wheel_next_res = len(coordinates)
            elif wheel_next_res < 1:
                wheel_next_res += len(coordinates)-1

            # line thickness
            thickness = 6*((i/len(coordinates)/1.2))
            thicknessNext = 6*(((i+1)/len(coordinates)/1.2))
            thicknessPrev = 6*(((i-1)/len(coordinates)/1.2))

            # find residue base points
            cur_points = self.ResiduePoints(cur_res, thickness, coordinates)
            next_points = self.ResiduePoints(next_res, thicknessNext, coordinates)
            prev_points = self.ResiduePoints(prev_res, thicknessPrev, coordinates)

            # lines
            cur_in = self.LineEquation(cur_points[1], cur_points[2])
            cur_out = self.LineEquation(cur_points[4], cur_points[3])
            next_in = self.LineEquation(next_points[1], next_points[2])
            next_out = self.LineEquation(next_points[4], next_points[3])
            prev_in = self.LineEquation(prev_points[1], prev_points[2])

            # line points
            p1 = cur_points[1];
            p2 = self.LineIntercept(cur_in['m'], cur_in['b'], next_out['m'], next_out['b'])
            p3 = self.LineIntercept(cur_in['m'], cur_in['b'], next_in['m'], next_in['b'])
            p4 = self.LineIntercept(cur_in['m'], cur_in['b'], prev_in['m'], prev_in['b'])
            p5 = self.LineIntercept(cur_out['m'], cur_out['b'], prev_in['m'], prev_in['b'])
            p6 = cur_points[4]

            # shorter line for the last residue
            if i == len(coordinates):
                move_in = self.MoveAlongLine(40,cur_in['m'],False,cur_in['x'],cur_in['y']);
                p4 = {'x':p1['x']+move_in['x'], 'y':p1['y']+move_in['y']}
                p5 = {'x':p6['x']+move_in['x'], 'y':p6['y']+move_in['y']}

            # draw line
            points = [p2,p1,p6,p5,p4,p3]
            points_txt = ""
            for coord in points:
                points_txt += str(coord['x'])+","+str(coord['y'])+" "

            # gradient
            gradientId = uniqid()
            angle = 1/tan(self.deg2rad(cur_in['m']))
            output +=  """
                <defs>
                    <linearGradient id='{}' x1='100%' y1='0%' x2='0%' y2='0%' gradientTransform='rotate({})'>
                    <stop offset='0%' stop-color='#00cc00' stop-opacity='1'/>
                    <stop offset='100%' stop-color='#006600' stop-opacity='1'/>
                    </linearGradient>
                </defs>
                """.format(gradientId,angle)

            lineFill = "white"

            # add SVG to output
            output += "<polyline stroke='black' stroke-width='0.5' points='"+points_txt+"' fill='"+lineFill+"'/>"

        return output;

    def LineIntercept(self,m1, b1, m2, b2):
        # line intercept
        intercept = {}
        intercept['x'] = (b2-b1)/(m1-m2)
        intercept['y'] = m1*intercept['x']+b1

        return intercept

    #Finds points on both sides of the residue center and the equation of a line to a reference residue
    def ResiduePoints(self, cur_res, thickness, coordinates):
        # FIXME write a general formula
        ref_residues = {
            20 : 2,
            16 : 2,
            12 : 2, # exception
            8 : 17,
            4 : 13,
            19 : 1,
            15 : 1,
            11 : 20,
            7 : 20,
            3 : 16,
            18 : 4,
            14 : 4,
            10 : 19,
            6 : 19,
            2 : 15,
            17 : 3,
            13 : 3,
            9 : 18,
            5 : 18,
            1 : 14}

        ori = {}

        ref_res = ref_residues[cur_res]

        # next-in-wheel residue orientation
        numSides = 4
        wheel_next_res = cur_res-numSides
        if wheel_next_res <= (1-numSides): # if current res is the last res, next is the first res
            wheel_next_res = len(coordinates)
        elif wheel_next_res < 1:
            wheel_next_res += len(coordinates)-1

        # point orientation
        if (coordinates[cur_res]['y'] >= coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] >
            coordinates[wheel_next_res]['x']):
            # SW
            ori['x'] = -1
            ori['y'] = -1

        elif (coordinates[cur_res]['y'] > coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] <
            coordinates[wheel_next_res]['x']):
            # NW
            ori['x'] = 1
            ori['y'] = 1

        elif (coordinates[cur_res]['y'] < coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] >
            coordinates[wheel_next_res]['x']):
            # SE
            ori['x'] = -1
            ori['y'] = -1

        elif (coordinates[cur_res]['y'] <= coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] <
            coordinates[wheel_next_res]['x']):
            # NE
            ori['x'] = 1
            ori['y'] = 1


        # slope of line from current residue to reference residue
        m = (coordinates[ref_res]['y']-coordinates[cur_res]['y'])/(coordinates[ref_res]['x']-
            coordinates[cur_res]['x'])

        # calculate coordinates of perpendicular line
        per_move = self.MoveAlongLine(thickness/2,m,True); # move thickness/2 pixels along a perpendicular line

        # define points
        points = [0];
        points.append({'x':coordinates[cur_res]['x']+per_move['x']*ori['x'], 'y':coordinates[cur_res]['y']+
            per_move['y']*ori['y']})
        points.append({'x':coordinates[ref_res]['x']+per_move['x']*ori['x'], 'y':coordinates[ref_res]['y']+
            per_move['y']*ori['y']})
        points.append({'x':points[2]['x']+per_move['x']*ori['x']*-2, 'y':points[2]['y']+per_move['y']*
            ori['y']*-2})
        points.append({'x':points[1]['x']+per_move['x']*ori['x']*-2, 'y':points[1]['y']+per_move['y']*
            ori['y']*-2})

        return points