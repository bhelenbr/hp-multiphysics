/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *  ShowMesh - SDL2 port
 *
 *  Cross-platform (macOS + Linux, and Windows if you add SDL2 there too)
 *  rewrite of the original X11/Xlib version of ShowMesh by Bojan Niceno.
 *  All the mesh-file parsing and mesh math is unchanged; every Xlib call
 *  has been replaced with SDL2 + SDL2_ttf equivalents.
 *
 *  Original author: Bojan NICENO (niceno@univ.trieste.it)
 *  SDL2 port: rewritten to drop the X11 dependency.
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#include <SDL.h>
#include <SDL_ttf.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdarg.h>

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
#define SMALL 1e-30
#define GREAT 1e+30

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Definitions for the mesh
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int Nn = 0, Ne = 0, Ns = 0;
double xmax =-GREAT, ymax=-GREAT, xmin=GREAT, ymin=GREAT, scl, X0, Y0;

struct ele
{
    int i,  j,  k;
    int mark;
    double xv, yv;
}
* elem;


struct sid
{
    int ea, eb;           /* left and right element */
    int a, b, c, d;       /* left, right, start and end point */
    int mark;             /* is it off, is on the boundary */
}
* side;


struct nod
{
    double x, y;
    int mark;
}
* node;
/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define MAIN_WDTH     750
#define MAIN_HGHT     500
#define BUTTON_HEIGHT  18
#define BUTTON_WIDTH   120
#define NBUTTONS       12
#define ON              0
#define OFF            -1

#define DELAUNAY  0
#define VORONOI   1
#define NODES     2
#define SIDES     3
#define ELEMENTS  4
#define MATERIALS 5
#define BOUNDARY  6
#define ZOOM      7
#define MOVE      8
#define FIT       9
#define QUIT      10
#define POSITION  11

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SDL / rendering state
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
static SDL_Window   *window   = NULL;
static SDL_Renderer *renderer = NULL;
static TTF_Font      *text_font = NULL;   /* button captions   */
static TTF_Font      *numb_font = NULL;   /* mesh numbering    */

static char *prog_name, *file_name;

static int main_wdth = MAIN_WDTH, main_hght = MAIN_HGHT;
static int draw_wdth, draw_hght;
static const int draw_x0 = 10, draw_y0 = 10; /* drawing area offset inside window */

static const SDL_Color COL_BLACK = {0, 0, 0, 255};
static const SDL_Color COL_WHITE = {255, 255, 255, 255};
static const SDL_Color COL_DRAG  = {200, 0, 0, 255};

struct butt_data
{
    int    x0, y0, hght, wdth, border;
    char   *caption;
    int    pressed;
}

butt_data[NBUTTONS]=
{
    { MAIN_WDTH-BUTTON_WIDTH-20,  10, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Delaunay",  ON},
    { MAIN_WDTH-BUTTON_WIDTH-20,  40, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Voronoi",   OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20,  70, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Nodes",     OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 100, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Sides",     OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 130, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Elements",  OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 160, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Materials", OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 190, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Boundary",  OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 220, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Zoom",      OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 250, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Move",      OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 280, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Fit",       OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 310, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Quit",      OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 340, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Position",  OFF},
};

static char position_buf[100] = "";

/*========================================================================*/
/*  Font discovery - no fonts are bundled, we look for common system      *
 *  fonts on macOS / Linux, or honour the SHOWMESH_FONT env var.          */
/*========================================================================*/
static const char *candidate_fonts[] =
{
    /* Linux (Debian/Ubuntu/Fedora/Arch common paths) */
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
    /* macOS */
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/Library/Fonts/Arial Bold.ttf",
    "/System/Library/Fonts/Supplemental/Helvetica.ttc",
    "/System/Library/Fonts/Helvetica.ttc",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    NULL
};

static char *find_font(void)
{
    const char *env = getenv("SHOWMESH_FONT");
    if(env && access(env, R_OK) == 0)
        return strdup(env);

    for(int i = 0; candidate_fonts[i] != NULL; i++)
        if(access(candidate_fonts[i], R_OK) == 0)
            return strdup(candidate_fonts[i]);

    return NULL;
}

/*========================================================================*/
static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    exit(-1);
}

/*========================================================================*/
static void load_fonts(void)
{
    char *font_path = find_font();

    if(font_path == NULL)
    {
        die("%s: could not find a TrueType font to use.\n"
            "  Set the SHOWMESH_FONT environment variable to a .ttf file, e.g.:\n"
            "    export SHOWMESH_FONT=/path/to/some/Bold.ttf\n"
            "  On Linux:  sudo apt install fonts-dejavu-core   (or fonts-liberation)\n"
            "  On macOS:  Arial Bold ships with the OS; if missing, point\n"
            "             SHOWMESH_FONT at any .ttf/.ttc font you have.\n",
            prog_name);
    }

    text_font = TTF_OpenFont(font_path, 13);
    numb_font = TTF_OpenFont(font_path, 11);

    if(text_font == NULL || numb_font == NULL)
        die("%s: cannot open font '%s': %s\n", prog_name, font_path, TTF_GetError());

    free(font_path);
}

/*========================================================================*/
/*  Small text helpers (replace XTextWidth / XDrawString)                 */
/*========================================================================*/
static void text_size(TTF_Font *font, const char *str, int *w, int *h)
{
    if(str == NULL || str[0] == '\0') { *w = 0; *h = TTF_FontHeight(font); return; }
    TTF_SizeText(font, str, w, h);
}

static void draw_text(TTF_Font *font, const char *str, int x, int y, SDL_Color color)
{
    if(str == NULL || str[0] == '\0') return;

    SDL_Surface *surf = TTF_RenderText_Blended(font, str, color);
    if(surf == NULL) return;

    SDL_Texture *tex = SDL_CreateTextureFromSurface(renderer, surf);
    SDL_Rect dst = { x, y, surf->w, surf->h };
    SDL_RenderCopy(renderer, tex, NULL, &dst);
    SDL_DestroyTexture(tex);
    SDL_FreeSurface(surf);
}

/*========================================================================*/
/*  Line drawing helpers (replace DrawLine / dashed GC / thick GC)        */
/*========================================================================*/
/* Bounding-box overlap test: is any part of segment (xl,yl)-(xr,yr),      *
 * given in mesh coordinates, inside the currently visible mesh-space     *
 * window? (The original X11 version tested each endpoint independently  *
 * against the same edge, which wrongly culled segments that poked past  *
 * *opposite* edges of the view - fixed here.)                            */
static int seg_visible(double xl, double yl, double xr, double yr)
{
    double x_lo = min(xl, xr), x_hi = max(xl, xr);
    double y_lo = min(yl, yr), y_hi = max(yl, yr);
    double vis_x_lo = -X0/scl,          vis_x_hi = (-X0 + draw_wdth)/scl;
    double vis_y_lo = (Y0 - draw_hght)/scl, vis_y_hi = Y0/scl;

    if(x_hi < vis_x_lo || x_lo > vis_x_hi) return 0;
    if(y_hi < vis_y_lo || y_lo > vis_y_hi) return 0;
    return 1;
}

/* mesh (x,y) -> window pixel coords (inside the drawing area) */
static void mesh_to_screen(double x, double y, int *sx, int *sy)
{
    *sx = (int)(x*scl + X0) + draw_x0;
    *sy = (int)(-y*scl + Y0) + draw_y0;
}

static void draw_mesh_line(double xl, double yl, double xr, double yr, SDL_Color color)
{
    int ixl, iyl, ixr, iyr;

    if(!seg_visible(xl, yl, xr, yr)) return;

    mesh_to_screen(xl, yl, &ixl, &iyl);
    mesh_to_screen(xr, yr, &ixr, &iyr);

    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
    SDL_RenderDrawLine(renderer, ixl, iyl, ixr, iyr);
}

static void draw_mesh_line_thick(double xl, double yl, double xr, double yr, SDL_Color color)
{
    int ixl, iyl, ixr, iyr;

    if(!seg_visible(xl, yl, xr, yr)) return;

    mesh_to_screen(xl, yl, &ixl, &iyl);
    mesh_to_screen(xr, yr, &ixr, &iyr);

    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
    SDL_RenderDrawLine(renderer, ixl,   iyl,   ixr,   iyr);
    SDL_RenderDrawLine(renderer, ixl+1, iyl,   ixr+1, iyr);
    SDL_RenderDrawLine(renderer, ixl-1, iyl,   ixr-1, iyr);
    SDL_RenderDrawLine(renderer, ixl,   iyl+1, ixr,   iyr+1);
    SDL_RenderDrawLine(renderer, ixl,   iyl-1, ixr,   iyr-1);
}

static void draw_mesh_line_dashed(double xl, double yl, double xr, double yr, SDL_Color color)
{
    int ixl, iyl, ixr, iyr;

    if(!seg_visible(xl, yl, xr, yr)) return;

    mesh_to_screen(xl, yl, &ixl, &iyl);
    mesh_to_screen(xr, yr, &ixr, &iyr);

    double dx = ixr - ixl, dy = iyr - iyl;
    double len = sqrt(dx*dx + dy*dy);
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
    if(len < 1.0) { SDL_RenderDrawPoint(renderer, ixl, iyl); return; }

    double ux = dx/len, uy = dy/len;
    const double dash = 4.0, gap = 4.0;
    double pos = 0.0;
    while(pos < len)
    {
        double seg_end = min(pos + dash, len);
        SDL_RenderDrawLine(renderer,
                            (int)(ixl + ux*pos),     (int)(iyl + uy*pos),
                            (int)(ixl + ux*seg_end), (int)(iyl + uy*seg_end));
        pos += dash + gap;
    }
}

/*========================================================================*/
static SDL_Rect button_rect(int b)
{
    SDL_Rect r;
    r.x = main_wdth - BUTTON_WIDTH - 20;
    r.y = butt_data[b].y0;
    r.w = butt_data[b].wdth;
    r.h = butt_data[b].hght;
    return r;
}

static int point_in_rect(int x, int y, SDL_Rect r)
{
    return (x >= r.x && x < r.x+r.w && y >= r.y && y < r.y+r.h);
}

static int point_in_draw_area(int x, int y)
{
    return (x >= draw_x0 && x < draw_x0+draw_wdth && y >= draw_y0 && y < draw_y0+draw_hght);
}

/*========================================================================*/
static void draw_button(int b)
{
    SDL_Rect r = button_rect(b);
    SDL_Color fg, bg;

    if(butt_data[b].pressed == ON) { bg = COL_BLACK; fg = COL_WHITE; }
    else                           { bg = COL_WHITE; fg = COL_BLACK; }

    SDL_SetRenderDrawColor(renderer, bg.r, bg.g, bg.b, 255);
    SDL_RenderFillRect(renderer, &r);
    SDL_SetRenderDrawColor(renderer, fg.r, fg.g, fg.b, 255);
    SDL_RenderDrawRect(renderer, &r);

    int tw, th;
    text_size(text_font, butt_data[b].caption, &tw, &th);
    draw_text(text_font, butt_data[b].caption, r.x + (r.w-tw)/2, r.y + (r.h-th)/2, fg);
}

static void draw_buttons(void)
{
    for(int b = 0; b < NBUTTONS; b++)
        draw_button(b);
}

static void draw_mesh(void); /* defined below; used by rebuild_cache() */

/*========================================================================*/
/*  Cached-frame rendering.                                                *
 *                                                                          *
 *  SDL renderers are swap-chain double/triple buffered: SDL_RenderPresent *
 *  flips to a *different* backbuffer each call. Drawing only a small      *
 *  partial update (e.g. just the Position readout) and presenting that    *
 *  makes each buffer in the chain hold a different, inconsistent frame -  *
 *  which is exactly what showed up as screen flashing while the mouse     *
 *  moved, and gets worse the bigger the mesh is.                          *
 *                                                                          *
 *  Fix: render the (expensive, thousands-of-lines) mesh + buttons into an *
 *  off-screen target texture only when the view actually changes (button  *
 *  toggle, resize, pan/zoom applied). Every call to present_frame() then  *
 *  just blits that texture (cheap, independent of mesh size) and draws    *
 *  the live bits - the position readout, and the drag rubber-band - fresh *
 *  on top, so every single present is a complete, self-consistent frame.  */
static SDL_Texture *cache_tex = NULL;
static int cache_w = -1, cache_h = -1;

static void ensure_cache(void)
{
    if(cache_tex == NULL || cache_w != main_wdth || cache_h != main_hght)
    {
        if(cache_tex) SDL_DestroyTexture(cache_tex);
        cache_tex = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888,
                                       SDL_TEXTUREACCESS_TARGET, main_wdth, main_hght);
        cache_w = main_wdth;
        cache_h = main_hght;
    }
}

/* Repaint the static part of the scene (border + mesh + buttons) into    *
 * the cache texture. Call whenever the view changes.                     */
static void rebuild_cache(void)
{
    ensure_cache();
    SDL_SetRenderTarget(renderer, cache_tex);

    SDL_SetRenderDrawColor(renderer, COL_WHITE.r, COL_WHITE.g, COL_WHITE.b, 255);
    SDL_RenderClear(renderer);

    SDL_Rect border = { draw_x0-3, draw_y0-3, draw_wdth+6, draw_hght+6 };
    SDL_SetRenderDrawColor(renderer, COL_BLACK.r, COL_BLACK.g, COL_BLACK.b, 255);
    SDL_RenderDrawRect(renderer, &border);
    SDL_Rect border2 = { draw_x0-2, draw_y0-2, draw_wdth+4, draw_hght+4 };
    SDL_RenderDrawRect(renderer, &border2);

    SDL_Rect clip = { draw_x0, draw_y0, draw_wdth, draw_hght };
    SDL_RenderSetClipRect(renderer, &clip);
    draw_mesh();
    SDL_RenderSetClipRect(renderer, NULL);

    draw_buttons();

    SDL_SetRenderTarget(renderer, NULL);
}

/* Cheap, always-full-frame present: blit the cache, then draw the live   *
 * drag rubber-band (if any) and the position readout fresh on top.       */
static void present_frame(int drag_mode, int drag_x0, int drag_y0, int drag_x, int drag_y)
{
    SDL_SetRenderTarget(renderer, NULL);
    SDL_RenderCopy(renderer, cache_tex, NULL, NULL);

    if(drag_mode == 1) /* DRAG_MOVE */
    {
        SDL_SetRenderDrawColor(renderer, COL_DRAG.r, COL_DRAG.g, COL_DRAG.b, 255);
        SDL_RenderDrawLine(renderer, draw_x0+drag_x0, draw_y0+drag_y0, draw_x0+drag_x, draw_y0+drag_y);
    }
    else if(drag_mode == 2) /* DRAG_ZOOM */
    {
        SDL_Rect rb = { draw_x0 + min(drag_x0,drag_x), draw_y0 + min(drag_y0,drag_y),
                         abs(drag_x-drag_x0), abs(drag_y-drag_y0) };
        SDL_SetRenderDrawColor(renderer, COL_DRAG.r, COL_DRAG.g, COL_DRAG.b, 255);
        SDL_RenderDrawRect(renderer, &rb);
    }

    SDL_Rect r = button_rect(POSITION);
    SDL_SetRenderDrawColor(renderer, COL_WHITE.r, COL_WHITE.g, COL_WHITE.b, 255);
    SDL_RenderFillRect(renderer, &r);
    SDL_SetRenderDrawColor(renderer, COL_BLACK.r, COL_BLACK.g, COL_BLACK.b, 255);
    SDL_RenderDrawRect(renderer, &r);
    int tw, th;
    text_size(text_font, position_buf, &tw, &th);
    draw_text(text_font, position_buf, r.x + (r.w-tw)/2, r.y + (r.h-th)/2, COL_BLACK);

    SDL_RenderPresent(renderer);
}

/*=========================================================================*/
static void draw_mesh(void)
{
    int    e, n, s, ea, eb;
    double x, y, xc, yc, xd, yd, x1, y1, x2, y2;
    char   numb[80];
    int    f_wdth;
    const int f_hght = 9;

    /***********************
     *  Draw Delaunay Mesh  *
     ***********************/
    if(butt_data[DELAUNAY].pressed == ON)
        for(s = 0; s < Ns; s++)
            if(side[s].mark != OFF)
            {
                xc = node[side[s].c].x; yc = node[side[s].c].y;
                xd = node[side[s].d].x; yd = node[side[s].d].y;
                draw_mesh_line(xc, yc, xd, yd, COL_BLACK);
            }

    /**********************
     *  Draw Voronoi Mesh  *
     **********************/
    if(butt_data[VORONOI].pressed == ON)
        for(s = 0; s < Ns; s++)
            if(side[s].mark != OFF)
            {
                if((ea = side[s].ea) != OFF)
                { x1 = elem[ea].xv; y1 = elem[ea].yv; }
                else
                { x1 = 0.5*(node[side[s].c].x+node[side[s].d].x);
                  y1 = 0.5*(node[side[s].c].y+node[side[s].d].y); }

                if((eb = side[s].eb) != OFF)
                { x2 = elem[eb].xv; y2 = elem[eb].yv; }
                else
                { x2 = 0.5*(node[side[s].c].x+node[side[s].d].x);
                  y2 = 0.5*(node[side[s].c].y+node[side[s].d].y); }

                draw_mesh_line_dashed(x1, y1, x2, y2, COL_BLACK);
            }

    for(s = 0; s < Ns; s++)
        if(side[s].mark > 0) /* side is on the boundary */
        {
            xc = node[side[s].c].x; yc = node[side[s].c].y;
            xd = node[side[s].d].x; yd = node[side[s].d].y;
            draw_mesh_line_thick(xc, yc, xd, yd, COL_BLACK);
        }

    if(butt_data[MATERIALS].pressed == ON || butt_data[ELEMENTS].pressed == ON)
        for(e = 0; e < Ne; e++)
        {
            x = (node[elem[e].i].x + node[elem[e].j].x + node[elem[e].k].x)/3.0;
            y = (node[elem[e].i].y + node[elem[e].j].y + node[elem[e].k].y)/3.0;

            if(butt_data[MATERIALS].pressed == ON) sprintf(numb, "%d", elem[e].mark);
            if(butt_data[ELEMENTS].pressed  == ON) sprintf(numb, "%d", e);

            f_wdth = 0; { int hh; text_size(numb_font, numb, &f_wdth, &hh); }
            int sx, sy; mesh_to_screen(x, y, &sx, &sy);
            draw_text(numb_font, numb, sx - f_wdth/2, sy - f_hght/2, COL_BLACK);
        }

    if(butt_data[NODES].pressed == ON)
        for(n = 0; n < Nn; n++)
        {
            x = node[n].x; y = node[n].y;
            sprintf(numb, "%d", n);
            f_wdth = 0; { int hh; text_size(numb_font, numb, &f_wdth, &hh); }
            int sx, sy; mesh_to_screen(x, y, &sx, &sy);
            draw_text(numb_font, numb, sx - f_wdth/2, sy - 5 - f_hght, COL_BLACK);
        }

    if(butt_data[SIDES].pressed == ON)
        for(s = 0; s < Ns; s++)
        {
            x = 0.75*node[side[s].c].x + 0.25*node[side[s].d].x;
            y = 0.75*node[side[s].c].y + 0.25*node[side[s].d].y;

            sprintf(numb, "%d", s);
            f_wdth = 0; { int hh; text_size(numb_font, numb, &f_wdth, &hh); }
            int sx, sy; mesh_to_screen(x, y, &sx, &sy);
            draw_text(numb_font, numb, sx - f_wdth/2, sy - f_hght/2, COL_BLACK);
        }

    /*****************************
     *  Draw Boundary Conditions  *
     *****************************/
    if(butt_data[BOUNDARY].pressed == ON)
    {
        for(s = 0; s < Ns; s++)
            if(side[s].mark > 0)
            {
                x = 0.5*(node[side[s].c].x + node[side[s].d].x);
                y = 0.5*(node[side[s].c].y + node[side[s].d].y);

                sprintf(numb, "%d", side[s].mark);
                f_wdth = 0; { int hh; text_size(numb_font, numb, &f_wdth, &hh); }
                int sx, sy; mesh_to_screen(x, y, &sx, &sy);
                draw_text(numb_font, numb, sx - f_wdth/2, sy - f_hght/2, COL_BLACK);
            }

        for(n = 0; n < Nn; n++)
            if(node[n].mark > 0)
            {
                x = node[n].x; y = node[n].y;
                sprintf(numb, "%d", node[n].mark);
                f_wdth = 0; { int hh; text_size(numb_font, numb, &f_wdth, &hh); }
                int sx, sy; mesh_to_screen(x, y, &sx, &sy);
                draw_text(numb_font, numb, sx - f_wdth/2, sy - f_hght/2, COL_BLACK);
            }
    }
}
/*-draw_mesh--------------------------------------------------------------*/


#ifdef EASYMESH
/*========================================================================*/
void load_mesh()
{
 int n, s, e, len;
 int d1, d2, d3, d4, d5, d6;
 char dummy[80];
 FILE *in;

 strcat(file_name, ".n");
 len=strlen(file_name);
 
/*--------+
|  Nodes  |
+--------*/
 in=fopen(file_name, "r");
 if(in==NULL) 
  {fprintf(stderr, "%s: cannot open file: %s\n\n", prog_name, file_name); 
   fflush(stdout);
   exit(-1);}
   
 fscanf(in, "%d", &Nn);
 node=(struct nod *) calloc(Nn, sizeof(struct nod)); 
 if(node==NULL) 
  {fprintf(stderr, "%s: cannot allocate enough memory\n\n", prog_name); 
   fflush(stdout);
   exit(-1);}

 xmax = -GREAT; xmin = GREAT;
 ymax = -GREAT; ymin = GREAT;
 for(n=0; n<Nn; n++)
  {
   fscanf(in, "%s %lf %lf %d", dummy, &node[n].x, &node[n].y, &node[n].mark);
   xmax=max(xmax, node[n].x); ymax=max(ymax, node[n].y);
   xmin=min(xmin, node[n].x); ymin=min(ymin, node[n].y);
  }
 fclose(in);

/*-----------+
|  Elements  |
+-----------*/
 file_name[len-1]='e';
 in=fopen(file_name, "r");
 if(in==NULL)
  {fprintf(stderr, "%s: cannot open file: %s\n\n", prog_name, file_name); 
   fflush(stdout);
   exit(-1);}

 fscanf(in, "%d", &Ne);

 elem=(struct ele *) calloc(Ne, sizeof(struct ele)); 
 if(elem==NULL) 
  {fprintf(stderr, "%s: cannot allocate enough memory\n\n", prog_name); 
   fflush(stdout);
   exit(-1);}

 for(e=0; e<Ne; e++)
  {
   fscanf(in, "%s %d %d %d %d %d %d %d %d %d %lf %lf %d", 
               dummy, &elem[e].i, &elem[e].j, &elem[e].k,
                      &d1, &d2, &d3, &d4, &d5, &d6,
                      &elem[e].xv, &elem[e].yv, &elem[e].mark);
  }
 fclose(in);

/*--------+
|  Sides  |
+--------*/
 file_name[len-1]='s';
 in=fopen(file_name, "r");
 if(in==NULL)
  {fprintf(stderr, "%s: cannot open file: %s\n\n", prog_name, file_name); 
   fflush(stdout);
   exit(-1);}

 fscanf(in, "%d", &Ns);

 side=(struct sid *) calloc(Ns, sizeof(struct sid)); 
 if(node==NULL) 
  {fprintf(stderr, "%s: cannot allocate enough memory\n\n", prog_name); 
   fflush(stdout);
   exit(-1);}

 for(s=0; s<Ns; s++)
  {
   fscanf(in, dummy);     
   fscanf(in, "%s %d %d %d %d %d", 
                      dummy, &side[s].c, &side[s].d, &side[s].ea, &side[s].eb,
                             &side[s].mark);
  }
 fclose(in);

}
/*------------------------------------------------------------------------*/
#else
/*========================================================================*/
void load_mesh(char *file_name, int bnum)
{
    int i, j, n, s, e, nsbd, type[100], num[100];
    int d1;
    FILE *in;

    in = fopen(file_name, "r");
    if(in == NULL)
        die("%s: cannot open file: %s\n\n", prog_name, file_name);

    int readN, readS, readE;
    fscanf(in, "%*[^:]:%d%*[^:]:%d%*[^:]:%d", &readN, &readS, &readE);
    Nn += readN;
    Ns += readS;
    Ne += readE;

    struct nod *tempnode = node;
    node = (struct nod *) calloc(Nn, sizeof(struct nod));
    if(node == NULL)
        die("%s: cannot allocate enough memory\n\n", prog_name);

    if(Nn > readN)
    {
        for(n = 0; n < Nn-readN; ++n)
            node[n] = tempnode[n];
        free(tempnode);
    }

    for(n = Nn-readN; n < Nn; n++)
    {
        fscanf(in, "%*[^:]:%lf%lf", &node[n].x, &node[n].y);
        node[n].mark = 0;
        xmax = max(xmax, node[n].x); ymax = max(ymax, node[n].y);
        xmin = min(xmin, node[n].x); ymin = min(ymin, node[n].y);
    }

    struct sid *tempside = side;
    side = (struct sid *) calloc(Ns, sizeof(struct sid));
    if(side == NULL)
        die("%s: cannot allocate enough memory\n\n", prog_name);

    if(Ns > readS)
    {
        for(n = 0; n < Ns-readS; ++n)
            side[n] = tempside[n];
        free(tempside);
    }

    for(s = Ns-readS; s < Ns; s++)
    {
        fscanf(in, "%*[^:]:%d%d", &side[s].c, &side[s].d);
        side[s].c += Nn-readN;
        side[s].d += Nn-readN;
        side[s].mark = 0;
        side[s].ea = 0;
        side[s].eb = 0;
    }

    struct ele *tempelem = elem;
    elem = (struct ele *) calloc(Ne, sizeof(struct ele));
    if(elem == NULL)
        die("%s: cannot allocate enough memory\n\n", prog_name);

    if(Ne > readE)
    {
        for(n = 0; n < Ne-readE; ++n)
            elem[n] = tempelem[n];
        free(tempelem);
    }

    for(e = Ne-readE; e < Ne; e++)
    {
        fscanf(in, "%*[^:]:%d%d%d", &elem[e].i, &elem[e].j, &elem[e].k);
        elem[e].i += Nn-readN;
        elem[e].j += Nn-readN;
        elem[e].k += Nn-readN;
        elem[e].xv = 0.0;
        elem[e].yv = 0.0;
        elem[e].mark = bnum;
    }

    /* SIDE BOUNDARY INFO HEADER */
    fscanf(in, "%*[^:]:%d", &nsbd);
    if(nsbd > 100) die("error: too many boundaries\n");

    for(i = 0; i < nsbd; ++i)
    {
        fscanf(in, "%*[^:]:%d%*[^:]:%d", &type[i], &num[i]);
        for(j = 0; j < num[i]; ++j)
        {
            fscanf(in, "%*[^:]:%d", &d1);
            side[d1+Ns-readS].mark = type[i];
        }
    }

    /* VERTEX BOUNDARY INFO HEADER */
    fscanf(in, "%*[^:]:%d", &nsbd);
    if(nsbd > 100) die("error: too many boundaries\n");

    for(i = 0; i < nsbd; ++i)
    {
        fscanf(in, "%*[^:]:%d", &type[i]);
        fscanf(in, "%*[^:]:%d", &d1);
        node[d1+Nn-readN].mark = type[i];
    }
    fclose(in);
}
/*------------------------------------------------------------------------*/
#endif

/*========================================================================*/
static void fit_view(void)
{
    scl = min((0.9*(double)draw_hght)/(ymax-ymin+SMALL),
              (0.9*(double)draw_wdth)/(xmax-xmin+SMALL));
    X0 = draw_wdth*0.05 - xmin*scl;
    Y0 = draw_hght - draw_hght*0.05 + ymin*scl;
}

/*========================================================================*/
static void init_window(int argc, char **argv)
{
    if(SDL_Init(SDL_INIT_VIDEO) != 0)
        die("%s: cannot initialize SDL: %s\n", prog_name, SDL_GetError());

    if(TTF_Init() != 0)
        die("%s: cannot initialize SDL_ttf: %s\n", prog_name, TTF_GetError());

    load_fonts();

    char title[512];
    snprintf(title, sizeof(title), "ShowMesh 1.0 - %s", argv[1]);

    window = SDL_CreateWindow(title,
                               SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                               MAIN_WDTH, MAIN_HGHT,
                               SDL_WINDOW_RESIZABLE);
    if(window == NULL)
        die("%s: cannot create window: %s\n", prog_name, SDL_GetError());

    SDL_SetWindowMinimumSize(window, BUTTON_WIDTH+200, NBUTTONS*(BUTTON_HEIGHT+8)+40);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if(renderer == NULL)
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    if(renderer == NULL)
        die("%s: cannot create renderer: %s\n", prog_name, SDL_GetError());

    main_wdth = MAIN_WDTH;
    main_hght = MAIN_HGHT;
    draw_wdth = main_wdth - BUTTON_WIDTH - 40;
    draw_hght = main_hght - 20;

    fit_view();
}

/*========================================================================*/
static void print_usage(void)
{
    printf("\n*********************************************************");
    printf("\n****************                        *****************");
    printf("\n****************   PROGRAM:  ShowMesh   *****************");
    printf("\n****************                        *****************");
    printf("\n****************      version  1.0      *****************");
    printf("\n****************                        *****************");
    printf("\n****************  Author: Bojan NICENO  *****************");
    printf("\n**************** niceno@univ.trieste.it *****************");
    printf("\n****************                        *****************");
    printf("\n*********************************************************");
    printf("\n\nUsage:  showmesh  <NAME>");
    printf("\n\nShowMesh uses input files created with EasyMesh.\n\n");
}

/*========================================================================*/
int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        print_usage();
        exit(-1);
    }

    prog_name = argv[0];
    file_name = argv[1];

#ifdef EASYMESH
	load_mesh();
#else 
    for(int i = 1; i < argc; ++i) {
        load_mesh(argv[i], i-1);
    }
#endif

    init_window(argc, argv);

    enum { DRAG_NONE = 0, DRAG_MOVE = 1, DRAG_ZOOM = 2 };
    int drag_mode = DRAG_NONE;
    int drag_x0 = 0, drag_y0 = 0, drag_x = 0, drag_y = 0;

    rebuild_cache();
    present_frame(drag_mode, drag_x0, drag_y0, drag_x, drag_y);

    int running = 1;
    SDL_Event ev;

    while(running && SDL_WaitEvent(&ev))
    {
        switch(ev.type)
        {
            case SDL_QUIT:
                running = 0;
                break;

            case SDL_WINDOWEVENT:
                if(ev.window.event == SDL_WINDOWEVENT_RESIZED ||
                   ev.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                {
                    main_wdth = ev.window.data1;
                    main_hght = ev.window.data2;
                    draw_wdth = main_wdth - BUTTON_WIDTH - 40;
                    draw_hght = main_hght - 20;
                    rebuild_cache();
                    present_frame(drag_mode, drag_x0, drag_y0, drag_x, drag_y);
                }
                else if(ev.window.event == SDL_WINDOWEVENT_EXPOSED)
                {
                    rebuild_cache();
                    present_frame(drag_mode, drag_x0, drag_y0, drag_x, drag_y);
                }
                break;

            case SDL_MOUSEBUTTONDOWN:
            {
                if(ev.button.button != SDL_BUTTON_LEFT) break;

                int mx = ev.button.x, my = ev.button.y;
                int b_hit = -1;
                for(int b = 0; b < NBUTTONS; b++)
                    if(point_in_rect(mx, my, button_rect(b))) { b_hit = b; break; }

                if(b_hit == QUIT)
                {
                    running = 0;
                    break;
                }

                if(b_hit >= 0)
                {
                    butt_data[b_hit].pressed = (butt_data[b_hit].pressed == ON) ? OFF : ON;

                    if(b_hit == MATERIALS && butt_data[MATERIALS].pressed == ON)
                        butt_data[ELEMENTS].pressed = OFF;
                    if(b_hit == ELEMENTS && butt_data[ELEMENTS].pressed == ON)
                        butt_data[MATERIALS].pressed = OFF;

                    if(butt_data[DELAUNAY].pressed == OFF)
                    {
                        butt_data[ELEMENTS].pressed  = OFF;
                        butt_data[MATERIALS].pressed = OFF;
                    }
                    if(butt_data[DELAUNAY].pressed == OFF && butt_data[VORONOI].pressed == OFF)
                        butt_data[NODES].pressed = OFF;

                    if(b_hit == FIT)
                    {
                        fit_view();
                        butt_data[FIT].pressed = OFF;
                    }

                    rebuild_cache();
                    present_frame(drag_mode, drag_x0, drag_y0, drag_x, drag_y);
                }
                else if(point_in_draw_area(mx, my) && drag_mode == DRAG_NONE)
                {
                    if(butt_data[MOVE].pressed == ON)
                    {
                        drag_mode = DRAG_MOVE;
                        drag_x0 = drag_x = mx - draw_x0;
                        drag_y0 = drag_y = my - draw_y0;
                    }
                    else if(butt_data[ZOOM].pressed == ON)
                    {
                        drag_mode = DRAG_ZOOM;
                        drag_x0 = drag_x = mx - draw_x0;
                        drag_y0 = drag_y = my - draw_y0;
                    }
                }
                break;
            }

            case SDL_MOUSEMOTION:
            {
                int mx = ev.motion.x, my = ev.motion.y;

                if(drag_mode != DRAG_NONE)
                {
                    drag_x = mx - draw_x0;
                    drag_y = my - draw_y0;
                    /* cache still shows the pre-drag view; only the rubber-band overlay moves */
                    present_frame(drag_mode, drag_x0, drag_y0, drag_x, drag_y);
                }
                else if(point_in_draw_area(mx, my))
                {
                    double lx = (mx - draw_x0 - X0)/scl;
                    double ly = (my - draw_y0 - Y0)/scl;
                    snprintf(position_buf, sizeof(position_buf), "%4.2e,%4.2e", lx, ly);
                    butt_data[POSITION].caption = position_buf;
                    /* mesh unchanged - cheap cache blit, no mesh redraw */
                    present_frame(DRAG_NONE, 0, 0, 0, 0);
                }
                break;
            }

            case SDL_MOUSEBUTTONUP:
            {
                if(ev.button.button != SDL_BUTTON_LEFT) break;

                if(drag_mode == DRAG_MOVE)
                {
                    X0 += (drag_x - drag_x0);
                    Y0 += (drag_y - drag_y0);
                    drag_mode = DRAG_NONE;
                    butt_data[MOVE].pressed = OFF;
                    rebuild_cache();
                    present_frame(drag_mode, 0, 0, 0, 0);
                }
                else if(drag_mode == DRAG_ZOOM)
                {
                    if(abs(drag_x-drag_x0) > 2 && abs(drag_y-drag_y0) > 2)
                    {
                        double x0_fiz = (min(drag_x0, drag_x) - X0)/scl;
                        double y0_fiz = (min(drag_y0, drag_y) - Y0)/scl;
                        double scl_new = (min((double)draw_wdth/abs(drag_x0-drag_x),
                                               (double)draw_hght/abs(drag_y0-drag_y))) * scl;
                        scl = scl_new;
                        X0 = -x0_fiz*scl;
                        Y0 = -y0_fiz*scl;
                    }
                    drag_mode = DRAG_NONE;
                    butt_data[ZOOM].pressed = OFF;
                    rebuild_cache();
                    present_frame(drag_mode, 0, 0, 0, 0);
                }
                break;
            }
        } /* end switch */
    } /* end while */

    if(cache_tex) SDL_DestroyTexture(cache_tex);
    TTF_CloseFont(text_font);
    TTF_CloseFont(numb_font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
