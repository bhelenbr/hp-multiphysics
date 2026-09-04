import logging
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.WARNING)

import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

import os
import sys
if '-c' in sys.argv:
    tp.session.connect(host='localhost',port=7600)

class myPlt:
  def __init__(plt_, axFntSz_, axTtlFntSz_, lgndFntSz_):
    plt_.axFntSz = axFntSz_
    plt_.axTtlFntSz = axTtlFntSz_
    plt_.lgndFntSz = lgndFntSz_

def loadMyData(No_blocks_, data_dir_, data_name_):
    for index in range(0, No_blocks_):
        file_name = f"{data_name_}_b{index}.dat"
        file_path = os.path.join(data_dir_, file_name)
        try:
            tp.data.load_tecplot(file_path,assign_strand_ids = False)
            print(f"Loaded file: {file_path}")
        except tp.exception.TecplotException as e:
            print(f"Error loading {file_path}: {e}")

def myAxesStyle(plot_,mPlt_, xtitle_, ytitle_, xmin_, xmax_, ymin_, ymax_):
    plot_.axes.preserve_scale=True # This allows getting the entire domain
    # Ticks:
    plot_.axes.x_axis.ticks.direction = TickDirection.Out
    plot_.axes.x_axis.ticks.auto_spacing= False
    plot_.axes.x_axis.ticks.spacing = 1
    plot_.axes.x_axis.ticks.length = 0.5
    plot_.axes.x_axis.ticks.minor_length = 0.3
    plot_.axes.x_axis.ticks.line_thickness = 0.2
    plot_.axes.x_axis.tick_labels.font.typeface = 'Times'
    plot_.axes.x_axis.tick_labels.font.size_units = Units.Point
    plot_.axes.x_axis.tick_labels.font.size = mPlt_.axFntSz

    plot_.axes.y_axis.ticks.direction = TickDirection.Out
    plot_.axes.y_axis.ticks.auto_spacing = False
    plot_.axes.y_axis.ticks.spacing = 1
    plot_.axes.y_axis.ticks.length = 0.5
    plot_.axes.y_axis.ticks.minor_length=0.3
    plot_.axes.y_axis.ticks.line_thickness = 0.2
    plot_.axes.y_axis.tick_labels.font.typeface = 'Times'
    plot_.axes.y_axis.tick_labels.font.size_units = Units.Point
    plot_.axes.y_axis.tick_labels.font.size = mPlt_.axFntSz

    # Titles:
    plot_.axes.x_axis.title.title_mode = AxisTitleMode.UseText
    plot_.axes.x_axis.title.text = xtitle_
    plot_.axes.x_axis.title.font.typeface = 'Times'
    plot_.axes.x_axis.title.font.italic = True
    plot_.axes.x_axis.title.font.size_units = Units.Point
    plot_.axes.x_axis.title.font.size = mPlt_.axTtlFntSz

    plot_.axes.y_axis.title.title_mode = AxisTitleMode.UseText
    plot_.axes.y_axis.title.text = ytitle_
    plot_.axes.y_axis.title.font.typeface = 'Times'
    plot_.axes.y_axis.title.font.italic = True
    plot_.axes.y_axis.title.font.size_units = Units.Point
    plot_.axes.y_axis.title.font.size = mPlt_.axTtlFntSz
    plot_.axes.y_axis.title.offset = 3.5

    # Line properties:
    plot_.axes.x_axis.line.line_thickness = 0.3
    plot_.axes.y_axis.line.line_thickness = 0.3

    # Viewport:
    plot_.axes.viewport.left = 7.5
    plot_.axes.viewport.right = 87.5
    plot_.axes.viewport.top = 10
    plot_.axes.viewport.bottom = 30

    # Bounds:
    plot_.axes.x_axis.min = xmin_
    plot_.axes.x_axis.max = xmax_

    plot_.axes.y_axis.min = ymin_
    plot_.axes.y_axis.max = ymax_ # need it to be slightly larger than the domain to show the upper edge of domain

def myLegendStyle(legend_ , mPlt_, legendTitle_):
    legend_.header.use_custom_text = True
    legend_.header.text_type = TextType.LaTeX
    legend_.header.custom_text = legendTitle_
    legend_.header.font.size_units = Units.Point
    legend_.header.font.size = mPlt_.lgndFntSz

    legend_.number_font.typeface = 'Times'
    legend_.number_font.size_units = Units.Point
    legend_.number_font.size = mPlt_.lgndFntSz

    legend_.auto_resize = True
    legend_.row_spacing = 1.3
    legend_.box.box_type = tp.constant.TextBox.None_


with tp.session.suspend(): # this will speed up the code
    tp.new_layout()

    ###################################
    ###### Velocity Contour Plots #####
    ###################################

    ###########
    # 3-Phase #
    ###########
    pageV = tp.active_page()
    pageV.name = 'V Contours'

    # Load Data #
    #############
    data_dir = '/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/BaselineV12GARad/Results_sx1500/'
    data_name = 'data468'
    No_blocks = 3
    frame = tp.active_frame()
    frame.name = '3 Phase'
    frame.show_border = False
    loadMyData(No_blocks, data_dir, data_name)

    dataset = frame.dataset
    plot = frame.plot(PlotType.Cartesian2D)
    plot.show_edge = True

    # Define Equations #
    ####################
    k_l = 64 # thermal conductivity of Si melt [W/m.K]
    rho = 2530 # density of Si [kg/m^3]
    c_p = 1000 # specific heat of Si [J/kg.K]
    alpha = k_l/rho/c_p # thermal diffusivity [m^2/s]
    d = 0.01 # length scale [m]

    equation='{0} = sqrt({1}**2 + {2}**2)*{3}/{4}'.format('{Vmag}', '{V3}', '{V4}',alpha,d)
    tp.data.operate.execute_equation(equation)

    # Zone Names #
    ##############
    frame.activate()
    dataset.zone(0).name = 'Liquid'
    dataset.zone(1).name = 'Solid'
    dataset.zone(2).name = 'Gas'

    # Adjust Various Axis properties #
    ##################################
    mPlt = myPlt(12,12,10)
    xtitle = 'x<sub>1</sub> (cm)'
    ytitle = 'x<sub>2</sub> (cm)'
    xmin = -8.4
    xmax = 5.2
    ymin = -1.3
    ymax =  2.01 # need it to be slightly larger than the domain to show the upper edge of domain
    myAxesStyle(plot, mPlt,xtitle,ytitle, xmin, xmax, ymin, ymax)

    # Contour Settings #
    ####################
    plot.show_contour = True

    plot.contour(0).variable = dataset.variable('Vmag')
    # plot.contour(0).colormap_name = 'Sequential - Yellow/Orange/Red'
    plot.contour(0).colormap_name = 'Small Rainbow'
    plot.contour(0).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(0).colormap_filter.continuous_max = 20
    levels0 = [x for x in range(0,21,5)]
    # plot.contour(0).colormap_filter.continuous_max = 70
    # levels0 = [x for x in range(0,71,10)]
    plot.contour(0).levels.reset_levels(levels0)

    legend0 = plot.contour(0).legend
    legend0Title = '$V\ \mathrm{(m/s)}$'
    myLegendStyle(legend0, mPlt,legend0Title)
    legend0.vertical = False
    legend0.position = (88.5, 57.5) # for 0 to 20
    # legend0.position = (47, 57.5) # for 0 to 70

    plot.contour(1).variable = dataset.variable('Vmag')
    plot.contour(1).colormap_name = 'Small Rainbow'
    plot.contour(1).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(1).colormap_filter.continuous_max = 0.1
    levels1 = [x*0.01 for x in range(0,11,2)]
    plot.contour(1).levels.reset_levels(levels1)

    legend1 = plot.contour(1).legend
    legend1Title = '$V\ \mathrm{(m/s)}$'
    myLegendStyle(legend1, mPlt,legend1Title)
    legend1.vertical = False
    legend1.position=(89, 15)

    plot.fieldmap(dataset.zone('Gas')).contour.flood_contour_group = plot.contour(0)
    plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(1)
    plot.fieldmap(dataset.zone('Solid')).contour.show = False
    plot.fieldmap(dataset.zone('Solid')).shade.show = False

    # Streamlines #
    ###############
    plot.vector.u_variable = dataset.variable('V3')
    plot.vector.v_variable = dataset.variable('V4')
    plot.show_streamtraces = True
    # plot.streamtraces.timing.reset_delta()
    plot.streamtraces.add_rake(start_position=(xmax,0),end_position=(.9*xmax,ymax),stream_type=Streamtrace.TwoDLine,num_seed_points=6)
    plot.streamtraces.add_rake(start_position=(-xmax,0),end_position=(-.9*xmax,ymax),stream_type=Streamtrace.TwoDLine,num_seed_points=6)
    plot.streamtraces.line_thickness = 0.03
    plot.streamtraces.arrowhead_size = 0.6
    plot.streamtraces.step_size = .25
    plot.streamtraces.arrowhead_spacing = 33
    plot.streamtraces.max_steps = 6500

    plot.streamtraces.set_termination_line([(xmin,0),(xmin, ymin), (xmax,ymin), (xmax,0)])
    plot.streamtraces.termination_line.active = True
    plot.streamtraces.termination_line.show = False
    plot.streamtraces.add_rake(start_position=(0.9*xmax,ymin),end_position=(.9*xmax,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=6)
    plot.streamtraces.add_rake(start_position=(2.1,ymin),end_position=(2.1,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=3)
    plot.streamtraces.add_rake(start_position=(0,ymin),end_position=(0,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=3)
    plot.streamtraces.add(seed_point=[-1.5, -0.85], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-6, -0.3], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-6, -0.06], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-6, -1], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-7.5, -0.9], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-7.5, -1.1], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[0.54, -0.063], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-0.42, -0.013], stream_type=Streamtrace.TwoDLine)

    #####################
    # 3-Phase Zoomed-In #
    #####################

    # Add a rectangle shoing the zoomed-in region:
    tp.macro.execute_command('''$!AttachGeom 
      GeomType = Rectangle
      AnchorPos
        {
        X = -0.3
        Y = -0.2
        }
      Color = Custom2
      LineThickness = 0.1
      RawData
    0.63 0.45''')


    frame = pageV.add_frame(position=(5.05,3.54),size=(2.15,1.48))
    frame.show_border=False
    frame.transparent=True
    frame.name = '3 Phase - Zoomed In'

    loadMyData(No_blocks, data_dir, data_name)

    dataset = frame.dataset
    plot = frame.plot(PlotType.Cartesian2D)
    plot.show_edge = True

    equation='{0} = sqrt({1}**2 + {2}**2)*{3}/{4}'.format('{Vmag}', '{V3}', '{V4}',alpha,d)
    tp.data.operate.execute_equation(equation)

    dataset.zone(0).name = 'Liquid'
    dataset.zone(1).name = 'Solid'
    dataset.zone(2).name = 'Gas'

    plot.axes.x_axis.show = False
    plot.axes.y_axis.show = False
    plot.axes.x_axis.min = -0.3
    plot.axes.x_axis.max = 0.3
    plot.axes.y_axis.min = -0.2
    plot.axes.y_axis.max=0.25


    # Contour Settings #
    ####################
    plot.show_contour = True

    plot.contour(0).variable = dataset.variable('Vmag')
    # plot.contour(0).colormap_name = 'Sequential - Yellow/Orange/Red'
    plot.contour(0).colormap_name = 'Small Rainbow'
    plot.contour(0).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(0).colormap_filter.continuous_max = 70
    levels0 = [x for x in range(0,71,10)]
    plot.contour(0).levels.reset_levels(levels0)

    plot.contour(0).legend.show = False

    plot.contour(1).variable = dataset.variable('Vmag')
    plot.contour(1).colormap_name = 'Small Rainbow'
    plot.contour(1).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(1).colormap_filter.continuous_max = 0.1
    levels1 = [x*0.01 for x in range(0,11,2)]
    plot.contour(1).levels.reset_levels(levels1)

    plot.contour(1).legend.show = False

    plot.fieldmap(dataset.zone('Gas')).contour.flood_contour_group = plot.contour(0)
    plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(1)
    plot.fieldmap(dataset.zone('Solid')).contour.show = False
    plot.fieldmap(dataset.zone('Solid')).shade.show = False

    plot.fieldmaps(0,1,2).edge.line_thickness = 0.5

    # Streamlines #
    ###############
    plot.vector.u_variable = dataset.variable('V3')
    plot.vector.v_variable = dataset.variable('V4')
    plot.show_streamtraces = True
    # plot.streamtraces.timing.reset_delta()
    # plot.streamtraces.add_rake(start_position=(-0.3,0.02),end_position=(0.3,0.02),stream_type=Streamtrace.TwoDLine,num_seed_points=4)
    plot.streamtraces.add(seed_point=[0.1, -0.08], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-0.1, -0.08], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[0.08, 0.00], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-0.08, 0.00], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[0.06, -0.17], stream_type=Streamtrace.TwoDLine)
    plot.streamtraces.add(seed_point=[-0.14, -0.1], stream_type=Streamtrace.TwoDLine)

    plot.streamtraces.line_thickness = 0.2
    plot.streamtraces.arrowhead_size = 3
    plot.streamtraces.step_size = .25
    plot.streamtraces.arrowhead_spacing = 25
    plot.streamtraces.max_steps = 2500



tp.export.save_tiff('/Users/ali/Documents/HRG/GasEffectPaper/ThreePhaseV.tiff', width = 5000, region=ExportRegion.AllFrames, supersample=3)


#################
# 2-Phase Cases #
#################
class myString:
  def __init__(frme, name, fldr, dta):
    frme.name = name 
    frme.fldr = fldr 
    frme.dta = dta

string1  = myString("2Phase_h","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_h/Results_sx1500/","data175")
string2  = myString("2Phase_hS","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_hNss/Results_sx1500/","data415")
string3  = myString("2Phase_hSp","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_hNssNp/Results_sx1500/","data415")
stringList = [string1, string2, string3]

for i in range(3):
    with tp.session.suspend(): # this will speed up the code
        frame = pageV.add_frame()
        frame.name = stringList[i].name
        frame.show_border = False
        data_dir = stringList[i].fldr 
        data_name = stringList[i].dta
        No_blocks = 2
        loadMyData(No_blocks, data_dir, data_name)
        dataset = frame.dataset
        plot = frame.plot(PlotType.Cartesian2D)
        plot.show_edge = True
        # Define Equations #
        ####################
        equation='{0} = sqrt({1}**2 + {2}**2)*{3}/{4}'.format('{Vmag}', '{V3}', '{V4}',alpha,d)
        tp.data.operate.execute_equation(equation)
        
        # Zone Names #
        ##############
        dataset.zone(0).name = 'Liquid'
        dataset.zone(1).name = 'Solid'

        # Axis Settings #
        #################
        ymax =  0.01 # need it to be slightly larger than the domain to show the upper edge of domain
        myAxesStyle(plot, mPlt,xtitle,ytitle, xmin, xmax, ymin, ymax)
        plot.axes.preserve_scale=True # This allows getting the entire domain

        # Contour Settings #
        ####################
        plot.show_contour = True

        plot.contour(1).variable = dataset.variable('Vmag')
        plot.contour(1).colormap_name = 'Small Rainbow'
        plot.contour(1).colormap_filter.distribution = ColorMapDistribution.Continuous
        if i == 0:
            plot.contour(1).colormap_filter.continuous_max = 0.0025
            levels1 = [x*0.0001 for x in range(0,26,5)]
        else:
            plot.contour(1).colormap_filter.continuous_max = 0.1
            levels1 = [x*0.01 for x in range(0,11,2)]
            
        plot.contour(1).levels.reset_levels(levels1)

        legend1 = plot.contour(1).legend
        myLegendStyle(legend1, mPlt,legend1Title)
        legend1.vertical = False
        if i == 0:
            legend1.auto_resize = False
            plot.contour(1).labels.step = 2
            legend1.position=(90, 26.5)
        else:
            legend1.position=(39.5, 26.5)

        
        
        plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(1)
        plot.fieldmap(dataset.zone('Solid')).contour.show = False
        plot.fieldmap(dataset.zone('Solid')).shade.show = False

        # Streamlines #
        ###############
        plot.vector.u_variable = dataset.variable('V3')
        plot.vector.v_variable = dataset.variable('V4')
        plot.show_streamtraces = True
        # plot.streamtraces.timing.reset_delta()
        plot.streamtraces.line_thickness = 0.03
        plot.streamtraces.arrowhead_size = 0.6
        plot.streamtraces.step_size = .25
        plot.streamtraces.arrowhead_spacing = 33
        plot.streamtraces.max_steps = 6500

        if i == 0:
            plot.streamtraces.add_rake(start_position=(xmin,ymin),end_position=(xmin,0),stream_type=Streamtrace.TwoDLine,num_seed_points=6)
        else:
            plot.streamtraces.set_termination_line([(xmin,0),(xmin, ymin), (xmax,ymin), (xmax,0)])
            plot.streamtraces.termination_line.active = True
            plot.streamtraces.termination_line.show = False
            plot.streamtraces.add_rake(start_position=(0.9*xmax,ymin),end_position=(.9*xmax,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=6)
            plot.streamtraces.add_rake(start_position=(2.1,ymin),end_position=(2.1,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=3)
            plot.streamtraces.add_rake(start_position=(0,ymin),end_position=(0,0.2*ymin),stream_type=Streamtrace.TwoDLine,num_seed_points=3)
            plot.streamtraces.add(seed_point=[-1.5, -0.85], stream_type=Streamtrace.TwoDLine)
            plot.streamtraces.add(seed_point=[-6, -0.3], stream_type=Streamtrace.TwoDLine)
            plot.streamtraces.add(seed_point=[-6, -0.06], stream_type=Streamtrace.TwoDLine)
            plot.streamtraces.add(seed_point=[-6, -1], stream_type=Streamtrace.TwoDLine)
            plot.streamtraces.add(seed_point=[-7.5, -0.9], stream_type=Streamtrace.TwoDLine)
            plot.streamtraces.add(seed_point=[-7.5, -1.1], stream_type=Streamtrace.TwoDLine)

    tp.macro.execute_command('$!RedrawAll')
    tp.export.save_tiff('/Users/ali/Documents/HRG/GasEffectPaper/{}V.tiff'.format(frame.name), width = 5000, region=ExportRegion.AllFrames, supersample=3)

######################################
###### Temperature Contour Plots #####
######################################

with tp.session.suspend(): # this will speed up the code
    pageT = tp.add_page()
    pageT.name = 'T Contours'

    ###########
    # 3-Phase #
    ###########

    # Load Data #
    #############
    data_dir = '/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/BaselineV12GARad/Results_sx1500/'
    data_name = 'data468'
    No_blocks = 3
    frame = tp.active_frame()
    frame.name = '3 Phase'
    frame.show_border = False
    loadMyData(No_blocks, data_dir, data_name)

    dataset = frame.dataset
    plot = frame.plot(PlotType.Cartesian2D)
    plot.show_edge = True

    # Zone Names #
    ##############
    frame.activate()
    dataset.zone(0).name = 'Liquid'
    dataset.zone(1).name = 'Solid'
    dataset.zone(2).name = 'Gas'

    # Define Equations #
    ####################
    Tm = 1685 # Si melting temperature [K]
    equation1 = '{0} = {1}*{2}'.format('{Tdim}', '{V5}',Tm)
    tp.data.operate.execute_equation(equation1, zones=[dataset.zone('Liquid'),dataset.zone('Gas')])
    equation2 = '{0} = {1}*{2}'.format('{Tdim}', '{V3}',Tm)
    tp.data.operate.execute_equation(equation2, zones=[dataset.zone('Solid')])

    # Axis Settings #
    #################
    ymax = 2.01
    xmin =  -5 # Nothing is going on for smaller x. So, I will not plot it.
    mPlt = myPlt(9,9,8)
    myAxesStyle(plot, mPlt,xtitle,ytitle, xmin, xmax, ymin, ymax)

    # Contour Settings #
    ####################
    plot.show_contour = True

    plot.contour(0).variable = dataset.variable('Tdim')
    plot.contour(0).colormap_name = 'Small Rainbow'
    plot.contour(0).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(0).colormap_filter.continuous_max = 1710
    plot.contour(0).colormap_filter.continuous_min = 1660
    levels0 = [x for x in range(1660,1711,10)]
    plot.contour(0).levels.reset_levels(levels0)
    plot.contour(0).labels.step = 1

    legend0 = plot.contour(0).legend
    legend0Title = '$T\ \mathrm{(K)}$'
    myLegendStyle(legend0, mPlt,legend0Title)
    legend0.vertical = False
    legend0.auto_resize = False
    legend0.position = (84,15)

    plot.contour(1).variable = dataset.variable('Tdim')
    plot.contour(1).colormap_name = 'Small Rainbow'
    plot.contour(1).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(1).colormap_filter.continuous_max = 1700
    plot.contour(1).colormap_filter.continuous_min = 300
    levels1 = [x for x in range(300,1701,100)]
    plot.contour(1).levels.reset_levels(levels1)
    plot.contour(1).labels.step = 4

    legend1 = plot.contour(1).legend
    legend1Title = '$T\ \mathrm{(K)}$'
    myLegendStyle(legend1, mPlt,legend1Title)
    legend1.vertical =False
    legend1.auto_resize = False
    legend1.position=(48, 58)

    plot.fieldmap(dataset.zone('Gas')).contour.flood_contour_group = plot.contour(1)
    plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(0)
    plot.fieldmap(dataset.zone('Solid')).contour.flood_contour_group = plot.contour(0)

    plot.fieldmaps(0,1,2).contour.contour_type = ContourType.Overlay
    plot.fieldmaps(0,1,2).contour.line_thickness = 0.02

    plot.fieldmaps(dataset.zone('Liquid'),dataset.zone('Solid')).contour.line_group = plot.contour(0)
    plot.fieldmaps(dataset.zone('Gas')).contour.line_group = plot.contour(1)


    #####################
    # 3-Phase Zoomed-In #
    #####################

    # Add a rectangle shoing the zoomed-in region:
    tp.macro.execute_command('''$!AttachGeom 
      GeomType = Rectangle
      AnchorPos
        {
        X = -0.3
        Y = -0.2
        }
      Color = Custom2
      LineThickness = 0.1
      RawData
    0.63 0.45''')


    frame = pageT.add_frame(position=(5.05,3.54),size=(2.15,1.48))
    frame.show_border=False
    frame.transparent=True
    frame.name = '3 Phase - Zoomed In'

    loadMyData(No_blocks, data_dir, data_name)

    dataset = frame.dataset
    plot = frame.plot(PlotType.Cartesian2D)
    plot.show_edge = True

    dataset.zone(0).name = 'Liquid'
    dataset.zone(1).name = 'Solid'
    dataset.zone(2).name = 'Gas'

    equation1 = '{0} = {1}*{2}'.format('{Tdim}', '{V5}',Tm)
    tp.data.operate.execute_equation(equation1, zones=[dataset.zone('Liquid'),dataset.zone('Gas')])
    equation2 = '{0} = {1}*{2}'.format('{Tdim}', '{V3}',Tm)
    tp.data.operate.execute_equation(equation2, zones=[dataset.zone('Solid')])

    plot.axes.x_axis.show = False
    plot.axes.y_axis.show = False
    plot.axes.x_axis.min = -0.3
    plot.axes.x_axis.max = 0.3
    plot.axes.y_axis.min = -0.2
    plot.axes.y_axis.max = 0.25


    # Contour Settings #
    ####################
    plot.show_contour = True

    plot.contour(0).variable = dataset.variable('Tdim')
    plot.contour(0).colormap_name = 'Small Rainbow'
    plot.contour(0).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(0).colormap_filter.continuous_max = 1710
    plot.contour(0).colormap_filter.continuous_min = 1660
    levels0 = [x for x in range(1660,1711,10)]
    plot.contour(0).levels.reset_levels(levels0)

    plot.contour(0).legend.show = False

    plot.contour(1).variable = dataset.variable('Tdim')
    plot.contour(1).colormap_name = 'Small Rainbow'
    plot.contour(1).colormap_filter.distribution = ColorMapDistribution.Continuous
    plot.contour(1).colormap_filter.continuous_max = 1700
    plot.contour(1).colormap_filter.continuous_min = 300
    levels1 = [x for x in range(300,1701,100)]
    plot.contour(1).levels.reset_levels(levels1)

    plot.contour(1).legend.show = False

    plot.fieldmap(dataset.zone('Gas')).contour.flood_contour_group = plot.contour(1)
    plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(0)
    plot.fieldmap(dataset.zone('Solid')).contour.flood_contour_group = plot.contour(0)

    plot.fieldmaps(0,1,2).contour.contour_type = ContourType.Overlay
    plot.fieldmaps(0,1,2).contour.line_thickness = 0.1

    plot.fieldmaps(dataset.zone('Liquid'),dataset.zone('Solid')).contour.line_group = plot.contour(0)
    plot.fieldmaps(dataset.zone('Gas')).contour.line_group = plot.contour(1)

    plot.fieldmaps(0,1,2).edge.line_thickness = 0.5


tp.export.save_tiff('/Users/ali/Documents/HRG/GasEffectPaper/ThreePhaseT.tiff', width = 5000, region=ExportRegion.AllFrames, supersample=3)


#################
# 2-Phase Cases #
#################
class myString:
  def __init__(frme, name, fldr, dta):
    frme.name = name 
    frme.fldr = fldr 
    frme.dta = dta

string1  = myString("2Phase_h","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_h/Results_sx1500/","data175")
string2  = myString("2Phase_hS","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_hNss/Results_sx1500/","data415")
string3  = myString("2Phase_hSp","/Users/ali/Codes/Testing/tri_hp/buoyancy/ThreePhasePaper/TwoPhaseV3_GARad_hNssNp/Results_sx1500/","data415")
stringList = [string1, string2, string3]

for i in range(3):
    with tp.session.suspend(): # this will speed up the code
        frame = pageT.add_frame()
        frame.name = stringList[i].name
        frame.show_border = False
        data_dir = stringList[i].fldr 
        data_name = stringList[i].dta
        No_blocks = 2
        loadMyData(No_blocks, data_dir, data_name)
        dataset = frame.dataset
        plot = frame.plot(PlotType.Cartesian2D)
        plot.show_edge = True

        # Zone Names #
        ##############
        dataset.zone(0).name = 'Liquid'
        dataset.zone(1).name = 'Solid'
       
        # Define Equations #
        ####################
        equation2 = '{0} = {1}*{2}'.format('{Tdim}', '{V3}',Tm)
        tp.data.operate.execute_equation(equation2,zones=[dataset.zone('Solid')])
        equation1 = '{0} = {1}*{2}'.format('{Tdim}', '{V5}',Tm)
        tp.data.operate.execute_equation(equation1,zones=[dataset.zone('Liquid')])

        # Axis Settings #
        #################
        ymax = 0.01
        myAxesStyle(plot, mPlt,xtitle,ytitle, xmin, xmax, ymin, ymax)
        
        # Contour Settings #
        ####################
        plot.show_contour = True

        plot.contour(0).variable = dataset.variable('Tdim')
        plot.contour(0).colormap_name = 'Small Rainbow'
        plot.contour(0).colormap_filter.distribution = ColorMapDistribution.Continuous
        plot.contour(0).colormap_filter.continuous_max = 1710
        plot.contour(0).colormap_filter.continuous_min = 1660
        levels0 = [x for x in range(1660,1711,10)]
        plot.contour(0).levels.reset_levels(levels0)
        plot.contour(0).labels.step = 1


        legend0 = plot.contour(0).legend
        legend0Title = '$T\ \mathrm{(K)}$'
        myLegendStyle(legend0, mPlt,legend0Title)
        legend0.vertical = False
        legend0.auto_resize = False
        legend0.position = (50, 26)

        plot.fieldmap(dataset.zone('Liquid')).contour.flood_contour_group = plot.contour(0)
        plot.fieldmap(dataset.zone('Solid')).contour.flood_contour_group = plot.contour(0)
        plot.fieldmaps(0,1).contour.contour_type = ContourType.Overlay
        plot.fieldmaps(0,1).contour.line_thickness = 0.02

    tp.macro.execute_command('$!RedrawAll')
    tp.export.save_tiff('/Users/ali/Documents/HRG/GasEffectPaper/{}T.tiff'.format(frame.name), width = 5000, region=ExportRegion.AllFrames, supersample=3)

