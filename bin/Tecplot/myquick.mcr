
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