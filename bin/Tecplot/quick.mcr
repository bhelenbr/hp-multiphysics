#!MC 1000

##################################################################
#                                                                #
#            Default tecplot.mcr file.                           #
#                                                                #
#  This file is processed automatically by tecplot on startup.   #
#  The macro functions defined here will appear in the quick     #
#  macro panel (from the "Tools" menu).                          #
#                                                                #
##################################################################

#
# Function 3D Rotation Animation  - Ardith Ketchum, Revised 11/2003
#
# This macro function will rotate a 3D plot, allowing the user to 
# choose the axis for rotation, how many degrees the plot will be 
# rotated with each step, how many rotation steps will occur, and
# whether the animation will be saved to a file.  When a file is 
# saved, it will be found in the Tecplot home directory with an 
# appropriate name.  For example, an AVI file showing rotation 
# around the x axis will be called xaxisrotation.avi. 
#
$!Macrofunction Name = "3D Rotation Animation"
                ShowInMacroPanel = True
  $!PromptForTextString |RotationAxis|
    Instructions = "Enter axis for rotation (X,Y,Z,PSI,...)."
  $!PromptForTextString |RotationAngle|
    Instructions = "Enter number of degrees for each rotation step."
  $!PromptForTextString |NumSteps|
    Instructions = "Enter the number of rotation steps."
  $!PromptForTextString |Animation|
    Instructions = "Enter 0 for No Animation File, 1 for an AVI file, 2 for an RM file, 3 for a Flash file."
  $!If |Animation| != 0
    $!Varset |format| = "AVI"
    $!Varset |Extension| = "AVI"
    $!If |Animation| == 2
      $!Varset |format| = "Rastermetafile"
      $!Varset |Extension| = "RM"
    $!Elseif |Animation| == 3
      $!Varset |format| = "Flash"
      $!Varset |Extension| = "swf"

    $!Endif
    $!EXPORTSETUP EXPORTFORMAT = |format|
    $!EXPORTSETUP IMAGEWIDTH = 546
    $!EXPORTSETUP EXPORTFNAME = "|RotationAxis|AxisRotation.|Extension|"
    $!EXPORTSTART 
  $!Endif
  $!Loop |NumSteps|
    $!ROTATE3DVIEW |RotationAxis|
      ANGLE = |RotationAngle|
      ROTATEORIGINLOCATION = DEFINEDORIGIN
    $!Redraw
    $!If |Animation| != 0
      $!EXPORTNEXTFRAME 
    $!Endif
    $!Endloop
  $!If |Animation| != 0
    $!EXPORTFINISH 
  $!Endif
  $!Pause "Animation is completed.  If the rotated image is off-center or off-screen, reset the center of rotation and animate again."
$!Endmacrofunction

#
# Macro Function Reset Center of Rotation
# This allows Tecplot to change the center of rotation based on 
# the minimum and maximum values of x, y and z to allow better
# 3D Rotation.
#
$!Macrofunction Name = "Reset Center of Rotation"
  ShowInMacroPanel = True
$!RESET3DORIGIN ORIGINRESETLOCATION = DATACENTER
$!View Datafit
$!Endmacrofunction


#
# This macro function is used by the Cascade and Tile
# frames macros. It will prompt the user for PaperWidth
# and PaperHeight dimensions.  If the values are invalid,
# it will require the user to retype the values.
#
# In version 10 and later, |PAPERWIDTH| and |PAPERHEIGHT|
# are intrinsic macro variables, so this function is not
# needed, but it remains here for backward compatability
#
$!MACROFUNCTION
  NAME = "GetPaperDim"
  SHOWINMACROPANEL = FALSE  #We don't want this to show in the quick macro panel
  #
  # If you always use one paper size, instead of prompting
  # for the paper size, you can just set it explicitly.  Example:
  #
  # $!VARSET |PAPERWIDTH| = 11
  # $!VARSET |PAPERHEIGHT| = 8.5
  #
  
  #
  # Get the paper width
  #
  $!VARSET |PAPERDIMOK| = 0
  $!WHILE |PAPERDIMOK| == 0
    $!PROMPTFORTEXTSTRING |PAPERWIDTH|
      INSTRUCTIONS = "Enter the paper width."
    $!VARSET |PAPERDIMOK| = 1
    $!IF |PAPERWIDTH| < 1
      $!VARSET |PAPERDIMOK| = 0
      $!PAUSE "Paper width must be greater than or equal to 1."
    $!ENDIF
  $!ENDWHILE
  
  #
  # Get the paper height
  #
  $!VARSET |PAPERDIMOK| = 0
  $!WHILE |PAPERDIMOK| == 0
    $!PROMPTFORTEXTSTRING |PAPERHEIGHT|
      INSTRUCTIONS = "Enter the paper height."
    $!VARSET |PAPERDIMOK| = 1
    $!IF |PAPERHEIGHT| < 1
      $!VARSET |PAPERDIMOK| = 0
      $!PAUSE "Paper height must be greater than or equal to 1."
    $!ENDIF
  $!ENDWHILE 
  $!REMOVEVAR |PAPERDIMOK| 
$!ENDMACROFUNCTION


#
# Function Cascade Frames - Scott Fowler 11/2001
#
# This macro function will move and resize all the frames
# the layout in a cascading fashion. If there are too many
# frames to fit in one cascade, multiple layers of cascading
# frames will be produced. Frame order is maintained.
#
$!MACROFUNCTION
  NAME = "Cascade Frames"

  $!IF |TECPLOTVERSION| < 100
    $!RUNMACROFUNCTION "GetPaperDim"
  $!ENDIF

  # Change the delta to change the distance 
  # in which the frame are cascaded
  $!VARSET |DELTA|  = 0.25

  # These two values define the margin between the paper
  # and the edge of the outer frames.
  $!VARSET |STARTX| = 0.15
  $!VARSET |STARTY| = 0.15

  # This calculates a width and height such that all frames are
  # the same size.
  $!VARSET |FRAMEDIMOK| = 0
  $!VARSET |NUMLAYERS| = 1
  $!VARSET |FRAMESPERLAYER| = |NUMFRAMES|
  $!WHILE |FRAMEDIMOK| == 0
    $!VARSET |FRAMEWIDTH|  = ((|PAPERWIDTH|)  - ((|FRAMESPERLAYER|-1)*|DELTA|) - (|STARTX|*2))
    $!VARSET |FRAMEHEIGHT| = ((|PAPERHEIGHT|) - ((|FRAMESPERLAYER|-1)*|DELTA|) - (|STARTY|*2))

    $!VARSET |FRAMEDIMOK| = 1
    $!IF |FRAMEWIDTH| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEHEIGHT| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEDIMOK| == 0
      $!VARSET |NUMLAYERS| += 1
      $!VARSET |FRAMESPERLAYER| = (ceil(|NUMFRAMES|/|NUMLAYERS|))
    $!ENDIF
  $!ENDWHILE
  #
  # Now, reposition and resize each frame.
  #
  $!DRAWGRAPHICS NO
  $!VARSET |NUMFRAMESDRAWN| = 0
  $!LOOP |NUMLAYERS|
    $!VARSET |XPOS| = |STARTX|
    $!VARSET |YPOS| = |STARTY|
    $!LOOP |FRAMESPERLAYER|
      $!IF |NUMFRAMESDRAWN| < |NUMFRAMES|
        $!FRAMECONTROL POP
          FRAME = 1
        $!IF |LOOP| != 1
          $!VARSET |XPOS| += |DELTA|
          $!VARSET |YPOS| += |DELTA|
        $!ENDIF
        $!FRAMELAYOUT XYPOS
          {
            X = |XPOS|
            Y = |YPOS|
          }
        $!FRAMELAYOUT WIDTH  = |FRAMEWIDTH|
        $!FRAMELAYOUT HEIGHT = |FRAMEHEIGHT|
        $!VARSET |NUMFRAMESDRAWN| += 1
      $!ENDIF
    $!ENDLOOP
  $!ENDLOOP
  $!DRAWGRAPHICS YES

  $!REMOVEVAR |DELTA|
  $!REMOVEVAR |STARTX|
  $!REMOVEVAR |STARTY|
  $!REMOVEVAR |FRAMEDIMOK|
  $!REMOVEVAR |NUMLAYERS|
  $!REMOVEVAR |FRAMESPERLAYER|
  $!REMOVEVAR |FRAMEWIDTH|
  $!REMOVEVAR |FRAMEHEIGHT|
  $!REMOVEVAR |NUMFRAMESDRAWN|
  $!REMOVEVAR |XPOS|
  $!REMOVEVAR |YPOS|

  $!IF |TECPLOTVERSION| < 100
    $!REMOVEVAR |PAPERWIDTH|
    $!REMOVEVAR |PAPERHEIGHT|
  $!ENDIF

$!ENDMACROFUNCTION


#
# Function Tile Frames - Scott Fowler 11/2001
#                        Revised 8/2002
#                        Revised 8/15/2003 - Prompts for number of frames across rather than calculating automatically
#
# This macro function will move and resize all the frames
# the layout in a tile fashion.
#
# If there are too many frames using the specified |MARGIN|,
# the |MARGIN| will be reduced until all frames fit on the
# paper.  |MARGIN| can be negative which will result in 
# overlapping frames. Frame order is maintained.
#
$!MACROFUNCTION
  NAME = "Tile Frames"
  SHOWINMACROPANEL = TRUE

  # Before Tecplot 10 there were no |PAPERWIDTH| and |PAPERHEIGHT|
  # intrinsics, so we must prompt for them.
  $!IF |TECPLOTVERSION| < 100
    $!RUNMACROFUNCTION "GetPaperDim"
  $!ENDIF
  
  # Change the margin to change the distance 
  # between frames and the margin to the edge of the paper.
  # This value may be modified automatically if there are
  # too many frames for the specified area.
  $!VARSET |MARGIN|  = 0.0

  # These two values define the margin from the edge of the paper
  # to the edge of the outer frames.
  $!VARSET |XPOS|   = 0.15
  $!VARSET |YPOS|   = 0.15

  #
  # Get number of frames across
  #
  $!VARSET |NUMFRAMESACROSS| = 0
  $!WHILE |NUMFRAMESACROSS| == 0
    $!PROMPTFORTEXTSTRING |NUMFRAMESACROSS|
      INSTRUCTIONS = "Enter number of frames across paper."
    $!IF |NUMFRAMESACROSS| < 1
      $!VARSET |NUMFRAMESACROSS| = 0
      $!PAUSE "Number of frames across paper must be greater than or equal to 1."
    $!ENDIF
  $!ENDWHILE


  $!VARSET |TMP| = (int(|NUMFRAMESACROSS|))
  # If NUMFRAMESACROSS is not an integer, "round" up
  $!IF |NUMFRAMESACROSS| > |TMP|
    $!VARSET |NUMFRAMESACROSS| = (|TMP|)
    $!IF |PAPERWIDTH| >= |PAPERHEIGHT|
      $!VARSET |NUMFRAMESACROSS| += 1
    $!ENDIF
  $!ENDIF

  $!VARSET |NUMFRAMESDOWN| = (|NUMFRAMES|/|NUMFRAMESACROSS|)
  $!VARSET |TMP| = (int(|NUMFRAMESDOWN|))
  # If NUMFRAMESDOWN is not an integer, "round" up
  $!IF |NUMFRAMESDOWN| > |TMP|
    $!VARSET |NUMFRAMESDOWN| = (|TMP|+1)
  $!ENDIF

  # This calculates a width and height such that all frames are
  # the same size, and fit with the specified margin.  

  $!VARSET |FRAMEDIMOK| = 0
  $!WHILE |FRAMEDIMOK| == 0
    $!VARSET |FRAMEWIDTH|  = ((|PAPERWIDTH|-(|XPOS|*2)-((|NUMFRAMESACROSS|-1)*|MARGIN|))/(|NUMFRAMESACROSS|))
    $!VARSET |FRAMEHEIGHT| = ((|PAPERHEIGHT|-(|YPOS|*2)-((|NUMFRAMESDOWN|-1)*|MARGIN|))/(|NUMFRAMESDOWN|))

    $!VARSET |FRAMEDIMOK| = 1
    $!IF |FRAMEWIDTH| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEHEIGHT| < 0.5
      $!VARSET |FRAMEDIMOK| = 0
    $!ENDIF
    $!IF |FRAMEDIMOK| == 0
      $!VARSET |MARGIN| -= 0.01
    $!ENDIF
  $!ENDWHILE
 
  $!DRAWGRAPHICS NO

  $!VARSET |FRAMESDRAWN| = 1
  $!VARSET |NUMDOWNDRAWN| = 1

  $!LOOP |NUMFRAMESDOWN|        
    $!LOOP |NUMFRAMESACROSS|
      # Make sure that we want to draw this frame
      $!IF |FRAMESDRAWN| <= |NUMFRAMES|        
        $!FRAMECONTROL POP
          FRAME = 1
        $!VARSET |FRAMEX| = ((|LOOP|-1)*(|FRAMEWIDTH|+|MARGIN|)+|XPOS|)
        $!VARSET |FRAMEY| = ((|NUMDOWNDRAWN|-1)*(|FRAMEHEIGHT|+|MARGIN|)+|YPOS|)
        $!FRAMELAYOUT XYPOS
          {
            X = |FRAMEX|
            Y = |FRAMEY|
          }
        $!FRAMELAYOUT WIDTH  = |FRAMEWIDTH|
        $!FRAMELAYOUT HEIGHT = |FRAMEHEIGHT|
        $!VARSET |FRAMESDRAWN| += 1
      $!ENDIF
    $!ENDLOOP    
    # Increment the number of down rows.
    $!VARSET |NUMDOWNDRAWN| += 1
  $!ENDLOOP
  $!DRAWGRAPHICS YES

  $!REMOVEVAR |MARGIN|
  $!REMOVEVAR |XPOS|
  $!REMOVEVAR |YPOS|
  $!REMOVEVAR |NUMFRAMESACROSS|
  $!REMOVEVAR |NUMFRAMESDOWN|
  $!REMOVEVAR |TMP|
  $!REMOVEVAR |FRAMEDIMOK|
  $!REMOVEVAR |FRAMEWIDTH|
  $!REMOVEVAR |FRAMEHEIGHT|
  $!REMOVEVAR |FRAMESDRAWN|
  $!REMOVEVAR |NUMDOWNDRAWN|
  $!REMOVEVAR |FRAMEX|
  $!REMOVEVAR |FRAMEY|

  $!IF |TECPLOTVERSION| < 100
    $!REMOVEVAR |PAPERWIDTH|
    $!REMOVEVAR |PAPERHEIGHT|
  $!ENDIF

$!ENDMACROFUNCTION

# TO ADD TO LAUNCH SCRIPT:
# -qm $HOME/bin/quick.mcr
# TECPHYFILE=$HOME/Library/tecplot.phy
# export TECPHYFILE

# TO ADD TO .tecplot.cfg
# #!MC 1000
# $!INTERFACE 
#   USEOFFSCREENBITMAP = NO
#   MINPIXELSFORDRAG = 4
#   UNIXHELPBROWSERCMD = 'OPEN_RESULT=`open "@" 2#&1`; if test "$OPEN_RESULT" != "" ; then MOZILLA_RESULT=`mozilla -remote "OpenURL(@)" 2#&1`; if test "$MOZILLA_RESULT" = "No running window found." ; then mozilla "@"; fi; fi&' 
#   BEEPONFRAMEINTERRUPT = NO
#   SHOWFRAMEBORDERSWHENOFF = NO
#   MAXCUSTOMCOLORSININTERFACE = 55
#   OPENGLCONFIG
#     {
#     SCREENRENDERING
#       {
#       DEPTHBUFFERSIZE = 32
#       }
#     IMAGERENDERING
#       {
#       DEPTHBUFFERSIZE = 32
#       }
#     }
# $!FRAMELAYOUT 
#   SHOWBORDER = NO
#   XYPOS
#     {
#     X = 1
#     Y = 0.25
#     }
#   WIDTH = 9
#   HEIGHT = 8
# $!GLOBALFRAME 
#   SNAPTOGRID = NO
#   SNAPTOPAPER = NO
#   FRAMEHEADERHEIGHT = 0.2
#   FRAMEHEADERFORMAT = "&(FrameName) <math>=</math> &(Date) <math>=</math> &(DataSetTitle)" 
# $!FILECONFIG 
#   TEMPFILEPATH = "/Users/helenbrk/Library" 


$!MacroFunction
  Name = "Export Tiff"
  ShowInMacroPanel = True
$!DropDialog QuickMacroPanel
$!GETUSERINPUT |MFBD|
  INSTRUCTIONS = "filepath:" 
$!GETUSERINPUT |begnum|
  INSTRUCTIONS = "begin:"
$!GETUSERINPUT |endnum|
  INSTRUCTIONS = "end:"
$!GETUSERINPUT |incre|
  INSTRUCTIONS = "increment:"  
$!GETUSERINPUT |prefix|
  INSTRUCTIONS = "prefix:" 
$!GETUSERINPUT |suffix|
  INSTRUCTIONS = "suffix:" 
$!VARSET |count| = |begnum| 
$!WHILE |count| <= |endnum|
$!IF |count| < 10
$!VarSet |LFDSFN1| = '|MFBD|/|prefix|00|count||suffix|'
$!ENDIF
$!IF |count| >= 10
$!IF |count| < 100
$!VarSet |LFDSFN1| = '|MFBD|/|prefix|0|count||suffix|'
$!ENDIF
$!ENDIF
$!IF |count| >= 100
$!VarSet |LFDSFN1| = '|MFBD|/|prefix||count||suffix|'
$!ENDIF
$!READDATASET '|LFDSFN1|' 
  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYPOSITION
##########
#TEMPORARY
##########
#$!ACTIVEFIELDZONES -= [1]
#$!VIEW FIT
#$!TWODAXIS XDETAIL{RANGEMAX = 2.0}
#$!VIEW TRANSLATE
#  X = 0
#  Y = 10
#$!VIEW TRANSLATE
#  X = 0
#  Y = 10
#$!ACTIVEFIELDZONES += [1]
###########
$!REDRAWALL 
$!EXPORTSETUP EXPORTFORMAT = TIFF
$!EXPORTSETUP PALETTE = COLOR
$!VARSET |out| = |count|
$!IF |out| < 10
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/00|out|.tiff'
$!ENDIF
$!IF |out| >= 10
$!IF |out| < 100
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/0|out|.tiff'
$!ENDIF
$!ENDIF
$!IF |out| >= 100
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/|out|.tiff'
$!ENDIF
$!EXPORT 
  APPEND = NO
$!REMOVEVAR |LFDSFN1|
$!VARSET |count| += |incre|
$!ENDWHILE
$!REMOVEVAR |MFBD|
$!REMOVEVAR |begnum|
$!REMOVEVAR |endnum|
$!REMOVEVAR |incre|
$!REMOVEVAR |count|
$!REMOVEVAR |out|
$!EndMacroFunction

$!MacroFunction
  Name = "Export avi"
  ShowInMacroPanel = True
$!DropDialog QuickMacroPanel
$!GETUSERINPUT |MFBD|
  INSTRUCTIONS = "filepath:" 
$!EXPORTSETUP
EXPORTFORMAT = AVI
EXPORTFNAME = '|MFBD|/movie.avi'
$!EXPORTSTART
$!GETUSERINPUT |begnum|
  INSTRUCTIONS = "begin:"
$!GETUSERINPUT |endnum|
  INSTRUCTIONS = "end:"
$!GETUSERINPUT |incre|
  INSTRUCTIONS = "increment:"  
$!GETUSERINPUT |prefix|
  INSTRUCTIONS = "prefix:" 
$!GETUSERINPUT |suffix|
  INSTRUCTIONS = "suffix:" 
$!VARSET |count| = |begnum| 
$!WHILE |count| <= |endnum|
$!IF |count| < 10
$!VarSet |LFDSFN1| = '|MFBD|/|prefix|00|count||suffix|'
$!ENDIF
$!IF |count| >= 10
$!IF |count| < 100
$!VarSet |LFDSFN1| = '|MFBD|/|prefix|0|count||suffix|'
$!ENDIF
$!ENDIF
$!IF |count| >= 100
$!VarSet |LFDSFN1| = '|MFBD|/|prefix||count||suffix|'
$!ENDIF
$!READDATASET '|LFDSFN1|' 
  READDATAOPTION = NEW
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYPOSITION
$!REDRAWALL 
$!EXPORTNEXTFRAME
$!REMOVEVAR |LFDSFN1|
$!VARSET |count| += |incre|
$!ENDWHILE
$!EXPORTFINISH
$!REMOVEVAR |MFBD|
$!REMOVEVAR |begnum|
$!REMOVEVAR |endnum|
$!REMOVEVAR |incre|
$!REMOVEVAR |count|
$!EndMacroFunction

$!MacroFunction
  Name = "Zone Movie with Fit"
  ShowInMacroPanel = True
$!GETUSERINPUT |begnum|
  INSTRUCTIONS = "begin:"
$!GETUSERINPUT |endnum|
  INSTRUCTIONS = "end:"
$!GETUSERINPUT |incre|
$!EXPORTSETUP
EXPORTFORMAT = AVI
EXPORTFNAME = '|MFBD|/movie.avi'
$!EXPORTSTART
$!VARSET |count| = |begnum|
$!WHILE |count| <= |endnum|
$!ACTIVEFIELDZONES = [|count|]
$!VIEW FIT
$!REDRAWALL 
$!EXPORTNEXTFRAME
$!VARSET |count| += |incre|
$!ENDWHILE
$!EXPORTFINISH
$!EndMacroFunction


$!MacroFunction
  Name = "Convergence Plot"
  ShowInMacroPanel = True
$!DELETELINEMAPS 
$!LOOP |NUMZONES|
$!CREATELINEMAP 
$!LINEMAP [|LOOP|]  NAME = '&ZN&'
$!LINEMAP [|LOOP|]  ASSIGN{XAXISVAR = 1}
$!LINEMAP [|LOOP|]  ASSIGN{YAXISVAR = 3}
$!LINEMAP [|LOOP|]  ASSIGN{ZONE = |LOOP|}
$!LINEMAP [|LOOP|]  INDICES{IJKLINES = I}
$!ACTIVELINEMAPS += [|LOOP|]
$!ENDLOOP
$!XYLINEAXIS YDETAIL 1 {COORDSCALE = LOG}
$!VIEW FIT

$!EndMacroFunction