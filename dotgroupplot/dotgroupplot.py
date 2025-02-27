import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.transforms as transforms
from matplotlib.transforms import Affine2D
import numpy as np
import math
from matplotlib.patheffects import RendererBase
import matplotlib.patheffects as PathEffects
import seaborn as sns


def Rotate(v, theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]]) @ np.array(v).T

def AreCirclesOverlap(x1, y1, radius1, x2, y2, radius2, gap = 0):
    if ((x1-x2)**2 + (y1-y2)**2 < (radius1 + radius2 + gap)**2):
        return True
    return False

def DrawPolygon(x, y, radius, n_edge, facecolor, ax):
    """
    Draw a shape with n_edges centered at x,y, with raidus (center to corner)
    n_edge <= 0 means circle
    """
    if (n_edge > 0):
        points = []
        for j in range(n_edge):
            cx = 0
            cy = radius
            cx, cy = Rotate([cx, cy], (j-0.5) * 2*math.pi/n_edge)
            points.append([cx+x, cy+y])
        ax.add_patch(patches.Polygon(points, facecolor=facecolor, edgecolor="white",
                                    linewidth=0.01*radius))
    else:
        ax.add_patch(patches.Circle([x,y], radius=radius, facecolor=facecolor,
                                   edgecolor="white", linewidth=0.01*radius))
    return ax

def GetLapGap(radius, n_edge):
    """
    The gap between two laps for the same-group plotting
    """
    return max([11 / 5 * radius / math.cos(math.pi / n_edge),#  extra 20% in the radius for buffer
                      6 / 5 * radius / math.sin(math.pi / n_edge)] # gap should be large enough so dots in the same lap does not overlap
              )

def EstimateGroupRectangularRegion(center_x, center_y, data, n_edge, radius):
    """
    Estimate the rectangular region 
    
    Return:
    -----
    (x_min, y_min), (x_max, y_max)
    """
    count = len(data)
    lap_gap = GetLapGap(raidus, n_edge)
    n_lap = math.ceil(math.sqrt(n_edge**2 + 8*(count-1)*n_edge)/(2*n_edge) - 0.5 )
    return (center_x - lap_gap * n_lap - radius, \
            center_y - lap_gap * n_lap - radius), \
         (center_x +lap_gap * n_lap + radius, \
            center_y + lap_gap * n_lap + radius)

def EstimateGroupCircularRegion(center_x, center_y, data, n_edge, radius):
    """
    Esimate the circular region
    
    Return:
    -----
    center_x, center_y, radius
    """
    count = len(data)
    if (count == 1):
        return center_x, center_y, radius
    lap_gap = GetLapGap(radius, n_edge)
    n_lap = math.ceil(math.sqrt(n_edge**2 + 8*(count-1)*n_edge)/(2*n_edge) - 0.5 )
    return center_x, center_y, n_lap * lap_gap + radius

def CirclePushOut(x, y, r, x1, y1, r1, gap):
    """
    Push x, y, r away from x1,y1,r1 to make sure they have a gap at least gap
    
    Return:
    -----
    new center x,y
    """
    if (AreCirclesOverlap(x, y, r, x1, y1, r1, gap) == False):
        return x, y
    
    d = math.sqrt((x-x1)**2 + (y-y1)**2)
    norm = math.sqrt(x**2 + y**2)
    norm1 = math.sqrt(x1**2 + y1**2)
    
    # normalized vector for (x,y)
    nx, ny = x/norm, y/norm 
    
    # project (x1,y1) to v-(x,y)
    costheta = (nx * x1 + ny *y1) / norm1
    px, py = nx * norm1 * costheta, ny * norm1 * costheta
    
    projectd = math.sqrt((px - x1)**2 + (py - y1)**2)
    # push away this distance. r+r1+gap is the ideal center distance
    pushd = math.sqrt((r + r1 + gap) ** 2 - projectd**2)
    
    return nx * (pushd + norm1 * costheta), ny * (pushd + norm1 * costheta)

# Draw one group
def DrawGroup(center_x, center_y, data, hue, hue_order, hue_map, hue_style, facecolor, n_edge, radius, layout, ax):
    """
    Draw the plot for one group represent by data
    
    Return:
    -----
    axes object
    radius of the group circular region
    """
    
    count = len(data)
    coord = []
    
    region_radius = 0
    if (layout == "lap"):
        lap = 0
        
        #center difference between two laps
        lap_gap = GetLapGap(radius, n_edge)
               
        i = 0 
        while (i < count):
            tmpx = 0
            tmpy = lap * lap_gap 
            if (i == 0):
                tmpx, tmpy = Rotate([tmpx, tmpy], i * 2*math.pi/n_edge)
                coord.append([tmpx, tmpy])
                i += 1
            else:
                k = 0
                for j in range(i, min(i+lap*n_edge, count), lap):
                    tmpx1, tmpy1 = Rotate([tmpx, tmpy], k * 2*math.pi/n_edge)
                    tmpx2, tmpy2 = Rotate([tmpx, tmpy], (k+1) * 2*math.pi/n_edge)
                    for l in range(j, min(j + lap, count)):
                        coef = (j + lap - l) / lap
                        coord.append([coef*tmpx1 + (1-coef)*tmpx2, 
                                      coef*tmpy1 + (1-coef)*tmpy2])
                    k += 1
                i += lap * n_edge
            lap += 1
        region_radius = lap * lap_gap + radius
    elif (layout == "layer"):
        return ax
    elif (layout == "spiral"):
        return ax
    
    colors = []
    if (hue != None):
        if (hue_style == "vertical_slice"):
            # Reorder the points from left to right, top to top
            coord = sorted(coord, key=lambda x:[round(x[0],4), -round(x[1],4)])
        elif (hue_style == "horizontal_slice"):
            coord = sorted(coord, key=lambda x:[-round(x[1],4), round(x[0],4)])
        elif (hue_style == "pie_slice"):
            coord = sorted(coord, key=lambda x:[round(-math.atan2(round(x[0],4), round(-x[1],4)),4), round(x[0]**2 + x[1]**2, 4)])
        else: # "lap"
            coord = coord
            
        hue_counts = {h:cnt for h,cnt in data[hue].value_counts().items()}
        used = 0
        for h in hue_order:
            if (h not in hue_counts):
                continue
            tmp = hue_counts[h]
            used += tmp
            colors += [hue_map[h]] * tmp
        colors += [facecolor] * (count - used)
    else:
        colors = [facecolor] * count 
    
    #print(x, y)
    for i, (xi, yi) in enumerate(coord):
        DrawPolygon(xi+center_x, yi+center_y, radius, 0, colors[i], ax)
    return ax, region_radius

def dotgroupplot(data, group, hue = None, hue_order = None, facecolor = 'grey', 
                  layout="lap", n_edge=6, radius=1, hue_style="vertical_slice", group_firstlap_n = 6, 
                  palette=None, ax=None):
        """
        Draw the honeycomb plot 
        
        Args:
        -----
        data: panda dataframe holding the raw data
        group: the column name in data for grouping the rows.
        hue: color of the dot based on the column in the dataframe. Use default facecolor if it is not specified in this dictionary
        hue_order: order to assign the color; otherwise the levels are inferred from the data objects.
        facecolor: color of the dot if hue is not specified or not in hue_order.
        layout: how to draw the dot within one group. 
                "spiral": a shape of spiral. TODO
                "lap": a dot in the center, and other dots around it as laps. 
        n_edge: number of edge for each group's shape, default 6 is a like heagonal shape
        radius: dot size
        hue_style: how to add colors witin on group.
                "vertical_slice": from left to right
                "horizontal_slice": from top to bottom
                "pie_slice": similar to a pie chart
                "lap": in laps
        group_firstlap_n: how many groups to put in the first lap (default: 6)
        palette: seaborn's palette
        ax: axes objective
        
        Return:
        -----
        axes object
        """
        
        ax = ax or plt.gca()
        
        groupby = data.groupby(group)
        group_num = len(groupby)
        
        hue_map = {}
        if (hue is not None):    
            hlist = list(data[hue].unique())
            if (hue_order is None):
                hue_order = hlist
            mypalette = sns.color_palette(palette, len(hue_order))
            for i, h in enumerate(hue_order):
                hue_map[h] = mypalette[i]
                
        #if (hue != None):
        #    hue_categories = 
        i = 0 
        j = 0 # track within lap coordinate
        lap = 0
        lap_start = 0 
        lap_full = 1
        group_gap = 4 * radius
        xlim = [0, 0]
        ylim = [0, 0]
        group_regions = [[]] # a list of list, the first list is for each cyle, 
                            #and a second layer  element is a triple, centerx,y and raidus 
                            # The values here is before global shifting
        
        inner_angle = 2 * math.pi / group_firstlap_n # the turning angle for first lap
        shift_x = 0 # shift everything
        shift_y = 0
        for g, subdf in sorted(data.groupby("CDR3"), key=lambda item: len(item[1]), reverse=True):
            # find the appropriate positions
            center_x = 0
            center_y = 0 
            x, y, r = EstimateGroupCircularRegion(0, 0, subdf, n_edge, radius)
            if (lap > 0):
                if (lap == 1):
                    px, py, pr = group_regions[lap - 1][0]
                    center_x, center_y = Rotate([0, pr + group_gap + r], inner_angle * (i - lap_start))                    
                elif (lap > 1):
                    # The position of this group is determined by two in the inner layer and the previous one
                    k = math.floor(((i - lap_start)//lap) * (lap - 1)) # the index in the previous lap
                    px, py, pr = group_regions[lap - 1][k]
                    if ((i - lap_start)%lap != 0) :
                        # turning around the group 
                        # The other case is directly pushed a way from previous-lap group
                        px2, py2, pr2 = group_regions[lap - 1][k+lap-1]
                        turning = ((i-lap_start)%lap) * inner_angle / lap
                        #if ((px**2+py**2) > (px2**2 + py2**2)):
                        #    px, py, pr = px2, py2, pr2
                        #    turning = -turning
                        px, py = Rotate([px, py], turning)
                        
                        pl = math.sqrt(px**2 + py**2)
                        pl2 = math.sqrt(px2**2 + py2**2)
                        # Put a rough position, should be some where based on the two 
                        # groups in the previous lap
                        coef = ((i-lap_start)%lap)/lap
                        center_x = px / pl * ((1-coef)*pl + coef*pl2 + pr + r + group_gap)
                        center_y = py / pl * ((1-coef)*pl + coef*pl2 + pr + r + group_gap)   
                    else: # directly push
                        pl = math.sqrt(px**2 + py**2)
                        # Put a rough position
                        center_x = px / pl * (pl + pr + r + group_gap)
                        center_y = py / pl * (pl + pr + r + group_gap) 
                        
                    # push away from other dots in the previous lap
                    for px, py, pr in group_regions[lap - 1]:
                        center_x, center_y = CirclePushOut(center_x ,center_y, r,
                                px, py, pr, group_gap)
                # The previous groups on the same lap
                for px, py, pr in group_regions[lap]:
                    #px, py, pr = group_regions[lap][-1]
                    center_x, center_y = CirclePushOut(center_x ,center_y, r,
                                                  px, py, pr, group_gap)
            #DEBUG: test angles
            #if (i>0):
            #    print(i, math.acos(center_y / math.sqrt(center_x**2+center_y**2))/(2*math.pi)*360)        
            ax, group_radius = DrawGroup(center_x + shift_x, center_y + shift_x, subdf, hue, hue_order, hue_map, hue_style, facecolor, n_edge, radius, layout, ax)
            group_regions[lap].append([center_x, center_y, group_radius])
            
            xlim[0] = min(xlim[0], center_x - group_radius)
            xlim[1] = max(xlim[1], center_x + group_radius)
            ylim[0] = min(ylim[0], center_y - group_radius)
            ylim[1] = max(ylim[1], center_y + group_radius)
            
            #print([center_x, center_y, group_radius])
            i += 1
            j += 1
            if (i == lap_full):
                group_regions[lap].append( group_regions[lap][0] )
                group_regions.append([])
                lap += 1
                j = 0
                lap_start = i
                lap_full += lap * group_firstlap_n
                
        # Create the legend
        if (hue != None):
            legend_elements = []
            for h in hue_order:
                legend_elements.append(lines.Line2D([0], [0], marker='o', linestyle="None", 
                                                   color=hue_map[h], label=h, markersize=7))
            ax.legend(handles=legend_elements)
        
        xrange = xlim[1] - xlim[0]
        yrange = ylim[1] - ylim[0]
        xlim[0] -= xrange * 0.025
        xlim[1] += xrange * 0.025
        ylim[0] -= yrange * 0.025
        ylim[1] += yrange * 0.025
        ax.set(xlim=(xlim[0] + shift_x, xlim[1] + shift_x),
              ylim=(ylim[0] + shift_y, ylim[1] + shift_y))

        return ax