#include "video.h"

// Assumes polygons have no more than four sides and are clipped a
// maximum of four times by frustum. Must be increased for more sides
// or more clip planes.
#define MAX_POLY_VERTS 8

#define MAX_SCREEN_HEIGHT  2048
#define MOVEMENT_SPEED     3.0
#define VMOVEMENT_SPEED    3.0
#define MAX_MOVEMENT_SPEED 30.0
#define PI                 3.141592
#define ROLL_SPEED         (PI / 20.0)
#define PITCH_SPEED        (PI / 20.0)
#define YAW_SPEED          (PI / 20.0)
#define MAX_COORD          0x4000
#define NUM_FRUSTUM_PLANES 4
#define CLIP_PLANE_EPSILON 0.0001
#define MAX_SPANS          10000
#define MAX_SURFS          1000
#define MAX_EDGES          5000


typedef struct {
    double v[3];
} point_t;

typedef struct {
    double x;
    double y;
} point2D_t;

typedef struct {
    int x;
    int y;
    int count;
    int color;
} span_t;

typedef struct {
    double distance;
    point_t normal;
} plane_t;

typedef struct {
    int color;
    int numverts;
    point_t verts[MAX_POLY_VERTS];
    plane_t plane;
} polygon_t;

typedef struct surf_s {
    struct surf_s *pnext, *pprev;
    int color;
    int visxstart;
    double zinv00, zinvstepx, zinvstepy;
    int state;
} surf_t;

typedef struct {
    int numverts;
    point2D_t verts[MAX_POLY_VERTS];
} polygon2D_t;

typedef struct convexobject_s {
    struct convexobject_s* pnext;
    point_t center;
    int numpolys;
    polygon_t* ppoly;
} convexobject_t;

typedef struct edge_s {
    int x;
    int xstep;
    int leading;
    surf_t* psurf;
    struct edge_s *pnext, *pprev;
    struct edge_s* pnextremove;
} edge_t;

static SDL_bool quit = SDL_FALSE;

char szAppName[] = "Clip";           // The name of this application
char szTitle[] = "3D clipping demo"; // The title bar text
static double roll, pitch, yaw;
static double currentspeed;
static point_t currentpos;
static double fieldofview, xcenter, ycenter;
static double xscreenscale, yscreenscale, maxscale;
static double maxscreenscaleinv;
static int numobjects;
static double speedscale = 1.0;
static plane_t frustumplanes[NUM_FRUSTUM_PLANES];

static double mroll[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
static double mpitch[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
static double myaw[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
static point_t vpn, vright, vup;
static point_t xaxis = {1, 0, 0};
static point_t zaxis = {0, 0, 1};

static polygon_t polys0[] = {
    {15,
     4,
     {{-10, 10, -10}, {10, 10, -10}, {10, -10, -10}, {-10, -10, -10}},
     {10, {0, 0, -1}}},
    {14,
     4,
     {{10, 10, -10}, {10, 10, 10}, {10, -10, 10}, {10, -10, -10}},
     {10, {1, 0, 0}}},
    {13,
     4,
     {{10, 10, 10}, {-10, 10, 10}, {-10, -10, 10}, {10, -10, 10}},
     {10, {0, 0, 1}}},
    {12,
     4,
     {{-10, 10, 10}, {-10, 10, -10}, {-10, -10, -10}, {-10, -10, 10}},
     {10, {-1, 0, 0}}},
    {11,
     4,
     {{-10, 10, -10}, {-10, 10, 10}, {10, 10, 10}, {10, 10, -10}},
     {10, {0, 1, 0}}},
    {10,
     4,
     {{-10, -10, -10}, {10, -10, -10}, {10, -10, 10}, {-10, -10, 10}},
     {10, {0, -1, 0}}},
};

static polygon_t polys1[] = {
    {1,
     4,
     {{-200, 0, -200}, {-200, 0, 200}, {200, 0, 200}, {200, 0, -200}},
     {0, {0, 1, 0}}},
};

static polygon_t polys2[] = {
    {6,
     4,
     {{0, 10, 0}, {20, 10, 0}, {10, 10, -10}, {0, 10, -10}},
     {10, {0, 1, 0}}},
    {6,
     4,
     {{-10, 10, 10}, {0, 10, 10}, {0, 10, 0}, {-10, 10, 0}},
     {10, {0, 1, 0}}},
    {6,
     4,
     {{0, 10, 0}, {0, 10, -10}, {-10, 10, -10}, {-10, 10, 0}},
     {10, {0, 1, 0}}},
    {5,
     4,
     {{0, -10, 0}, {0, -10, -10}, {10, -10, -10}, {20, -10, 0}},
     {10, {0, -1, 0}}},
    {5,
     4,
     {{-10, -10, 10}, {-10, -10, 0}, {0, -10, 0}, {0, -10, 10}},
     {10, {0, -1, 0}}},
    {5,
     4,
     {{-10, -10, 0}, {-10, -10, -10}, {0, -10, -10}, {0, -10, 0}},
     {10, {0, -1, 0}}},
    {4,
     4,
     {{-10, 10, -10}, {10, 10, -10}, {10, -10, -10}, {-10, -10, -10}},
     {10, {0, 0, -1}}},
    {3,
     4,
     {{10, 10, -10}, {20, 10, 0}, {20, -10, 0}, {10, -10, -10}},
     {14.14, {0.707, 0, -0.707}}},
    {2,
     4,
     {{20, 10, 0}, {0, 10, 0}, {0, -10, 0}, {20, -10, 0}},
     {0, {0, 0, 1}}},
    {9,
     4,
     {{0, 10, 0}, {0, 10, 10}, {0, -10, 10}, {0, -10, 0}},
     {0, {1, 0, 0}}},
    {15,
     4,
     {{0, 10, 10}, {-10, 10, 10}, {-10, -10, 10}, {0, -10, 10}},
     {10, {0, 0, 1}}},
    {14,
     4,
     {{-10, 10, 10}, {-10, 10, -10}, {-10, -10, -10}, {-10, -10, 10}},
     {10, {-1, 0, 0}}},
};

extern convexobject_t objecthead;

static convexobject_t objects[] = {
    {&objects[1], {-50, 0, 70}, SDL_arraysize(polys2), polys2},
    {&objects[2], {0, 20, 70}, SDL_arraysize(polys0), polys0},
    {&objects[3], {50, 0, 70}, SDL_arraysize(polys0), polys0},
    {&objects[4], {-50, 0, -70}, SDL_arraysize(polys2), polys2},
    {&objects[5], {0, 20, -70}, SDL_arraysize(polys2), polys2},
    {&objects[6], {50, 30, -70}, SDL_arraysize(polys0), polys0},
    {&objects[7], {-50, 15, 0}, SDL_arraysize(polys0), polys0},
    {&objects[8], {50, 15, 0}, SDL_arraysize(polys2), polys2},
    {&objects[9], {0, 50, 0}, SDL_arraysize(polys2), polys2},
    {&objects[10], {-100, 100, 115}, SDL_arraysize(polys2), polys2},
    {&objects[11], {-100, 150, 120}, SDL_arraysize(polys0), polys0},
    {&objects[12], {100, 200, 100}, SDL_arraysize(polys0), polys0},
    {&objects[13], {100, 100, 100}, SDL_arraysize(polys2), polys2},
    {&objecthead, {0, -20, 0}, SDL_arraysize(polys1), polys1},
};

// Head and tail for the object list
convexobject_t objecthead = {&objects[0]};

// Span, edge, and surface lists
static span_t spans[MAX_SPANS];
static edge_t edges[MAX_EDGES];
static surf_t surfs[MAX_SURFS];

// Bucket list of new edges to add on each scan line
static edge_t newedges[MAX_SCREEN_HEIGHT];

// Bucket list of edges to remove on each scan line
static edge_t* removeedges[MAX_SCREEN_HEIGHT];

// Head and tail for the active edge list
static edge_t edgehead;
static edge_t edgetail;

// Edge used as sentinel of new edge lists
static edge_t maxedge = {0x7FFFFFFF};

// Head/tail/sentinel/background surface of active surface stack
static surf_t surfstack;

// pointers to next available surface and edge
static surf_t* pavailsurf;
static edge_t* pavailedge;

static int currentcolor;

static void UpdateWorld(void);

/////////////////////////////////////////////////////////////////////
// 3-D dot product.
/////////////////////////////////////////////////////////////////////
static double DotProduct(const point_t* vec1, const point_t* vec2) {
    return vec1->v[0] * vec2->v[0] + vec1->v[1] * vec2->v[1]
           + vec1->v[2] * vec2->v[2];
}

/////////////////////////////////////////////////////////////////////
// Concatenate two 3x3 matrices.
/////////////////////////////////////////////////////////////////////
static void MConcat(double in1[3][3], double in2[3][3], double out[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = in1[i][0] * in2[0][j] + in1[i][1] * in2[1][j]
                        + in1[i][2] * in2[2][j];
        }
    }
}

/////////////////////////////////////////////////////////////////////
// Project viewspace polygon vertices into screen coordinates.
// Note that the y axis goes up in worldspace and viewspace, but
// goes down in screenspace.
/////////////////////////////////////////////////////////////////////
static void ProjectPolygon(const polygon_t* ppoly, polygon2D_t* ppoly2D) {
    for (int i = 0; i < ppoly->numverts; i++) {
        double zrecip = 1.0 / ppoly->verts[i].v[2];
        ppoly2D->verts[i].x =
            ppoly->verts[i].v[0] * zrecip * maxscale + xcenter;
        ppoly2D->verts[i].y =
            ycenter - (ppoly->verts[i].v[1] * zrecip * maxscale);
    }
    ppoly2D->numverts = ppoly->numverts;
}

/////////////////////////////////////////////////////////////////////
// Move the view position and set the world->view transform.
/////////////////////////////////////////////////////////////////////
static void UpdateViewPos(void) {
    int i;
    point_t motionvec;
    double s, c, mtemp1[3][3], mtemp2[3][3];

    // Move in the view direction, across the x-y plane, as if
    // walking. This approach moves slower when looking up or
    // down at more of an angle
    motionvec.v[0] = DotProduct(&vpn, &xaxis);
    motionvec.v[1] = 0.0;
    motionvec.v[2] = DotProduct(&vpn, &zaxis);

    for (i = 0; i < 3; i++) {
        currentpos.v[i] += motionvec.v[i] * currentspeed;
        currentpos.v[i] = SDL_clamp(currentpos.v[i], -MAX_COORD, MAX_COORD);
    }

    // Set up the world-to-view rotation.
    // Note: much of the work done in concatenating these matrices
    // can be factored out, since it contributes nothing to the
    // final result; multiply the three matrices together on paper
    // to generate a minimum equation for each of the 9 final elements
    s = sin(roll);
    c = cos(roll);
    mroll[0][0] = c;
    mroll[0][1] = s;
    mroll[1][0] = -s;
    mroll[1][1] = c;

    s = sin(pitch);
    c = cos(pitch);
    mpitch[1][1] = c;
    mpitch[1][2] = s;
    mpitch[2][1] = -s;
    mpitch[2][2] = c;

    s = sin(yaw);
    c = cos(yaw);
    myaw[0][0] = c;
    myaw[0][2] = -s;
    myaw[2][0] = s;
    myaw[2][2] = c;

    MConcat(mroll, myaw, mtemp1);
    MConcat(mpitch, mtemp1, mtemp2);

    // Break out the rotation matrix into vright, vup, and vpn.
    // We could work directly with the matrix; breaking it out
    // into three vectors is just to make things clearer
    for (i = 0; i < 3; i++) {
        vright.v[i] = mtemp2[0][i];
        vup.v[i] = mtemp2[1][i];
        vpn.v[i] = mtemp2[2][i];
    }

    // Simulate crude friction
    if (currentspeed > (MOVEMENT_SPEED * speedscale / 2.0))
        currentspeed -= MOVEMENT_SPEED * speedscale / 2.0;
    else if (currentspeed < -(MOVEMENT_SPEED * speedscale / 2.0))
        currentspeed += MOVEMENT_SPEED * speedscale / 2.0;
    else
        currentspeed = 0.0;
}

/////////////////////////////////////////////////////////////////////
// Rotate a vector from viewspace to worldspace.
/////////////////////////////////////////////////////////////////////
static void BackRotateVector(const point_t* pin, point_t* pout) {
    // Rotate into the world orientation
    for (int i = 0; i < 3; i++) {
        pout->v[i] = pin->v[0] * vright.v[i] + pin->v[1] * vup.v[i]
                     + pin->v[2] * vpn.v[i];
    }
}

/////////////////////////////////////////////////////////////////////
// Transform a point from worldspace to viewspace.
/////////////////////////////////////////////////////////////////////
static void TransformPoint(const point_t* pin, point_t* pout) {
    point_t tvert;

    // Translate into a viewpoint-relative coordinate
    for (int i = 0; i < 3; i++) {
        tvert.v[i] = pin->v[i] - currentpos.v[i];
    }

    // Rotate into the view orientation
    pout->v[0] = DotProduct(&tvert, &vright);
    pout->v[1] = DotProduct(&tvert, &vup);
    pout->v[2] = DotProduct(&tvert, &vpn);
}

/////////////////////////////////////////////////////////////////////
// Transform a polygon from worldspace to viewspace.
/////////////////////////////////////////////////////////////////////
static void TransformPolygon(const polygon_t* pinpoly, polygon_t* poutpoly) {
    for (int i = 0; i < pinpoly->numverts; i++) {
        TransformPoint(&pinpoly->verts[i], &poutpoly->verts[i]);
    }
    poutpoly->numverts = pinpoly->numverts;
}

/////////////////////////////////////////////////////////////////////
// Returns true if polygon faces the viewpoint, assuming a clockwise
// winding of vertices as seen from the front.
/////////////////////////////////////////////////////////////////////
static int PolyFacesViewer(const polygon_t* ppoly, const plane_t* pplane) {
    point_t viewvec;

    for (int i = 0; i < 3; i++) {
        viewvec.v[i] = ppoly->verts[0].v[i] - currentpos.v[i];
    }

    // Use an epsilon here so we don't get polygons tilted so
    // sharply that the gradients are unusable or invalid
    if (DotProduct(&viewvec, &pplane->normal) < -0.01) {
        return 1;
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////
// Set up a clip plane with the specified normal.
/////////////////////////////////////////////////////////////////////
static void SetWorldspaceClipPlane(point_t* normal, plane_t* plane) {
    // Rotate the plane normal into worldspace
    BackRotateVector(normal, &plane->normal);
    plane->distance =
        DotProduct(&currentpos, &plane->normal) + CLIP_PLANE_EPSILON;
}

/////////////////////////////////////////////////////////////////////
// Set up the planes of the frustum, in worldspace coordinates.
/////////////////////////////////////////////////////////////////////
static void SetUpFrustum(void) {
    double angle, s, c;
    point_t normal;

    angle = atan(2.0 / fieldofview * maxscale / xscreenscale);
    s = sin(angle);
    c = cos(angle);

    // Left clip plane
    normal.v[0] = s;
    normal.v[1] = 0;
    normal.v[2] = c;
    SetWorldspaceClipPlane(&normal, &frustumplanes[0]);

    // Right clip plane
    normal.v[0] = -s;
    SetWorldspaceClipPlane(&normal, &frustumplanes[1]);

    angle = atan(2.0 / fieldofview * maxscale / yscreenscale);
    s = sin(angle);
    c = cos(angle);

    // Bottom clip plane
    normal.v[0] = 0;
    normal.v[1] = s;
    normal.v[2] = c;
    SetWorldspaceClipPlane(&normal, &frustumplanes[2]);

    // Top clip plane
    normal.v[1] = -s;
    SetWorldspaceClipPlane(&normal, &frustumplanes[3]);
}

/////////////////////////////////////////////////////////////////////
// Clip a polygon to a plane.
/////////////////////////////////////////////////////////////////////
static int ClipToPlane(polygon_t* pin, const plane_t* pplane, polygon_t* pout) {
    int i, j, nextvert, curin, nextin;
    double curdot, nextdot, scale;
    point_t *pinvert, *poutvert;

    pinvert = pin->verts;
    poutvert = pout->verts;

    curdot = DotProduct(pinvert, &pplane->normal);
    curin = (curdot >= pplane->distance);

    for (i = 0; i < pin->numverts; i++) {
        nextvert = (i + 1) % pin->numverts;

        // Keep the current vertex if it's inside the plane
        if (curin)
            *poutvert++ = *pinvert;

        nextdot = DotProduct(&pin->verts[nextvert], &pplane->normal);
        nextin = (nextdot >= pplane->distance);

        // Add a clipped vertex if one end of the current edge is
        // inside the plane and the other is outside
        if (curin != nextin) {
            scale = (pplane->distance - curdot) / (nextdot - curdot);
            for (j = 0; j < 3; j++) {
                poutvert->v[j] =
                    pinvert->v[j]
                    + ((pin->verts[nextvert].v[j] - pinvert->v[j]) * scale);
            }
            poutvert++;
        }

        curdot = nextdot;
        curin = nextin;
        pinvert++;
    }

    pout->numverts = poutvert - pout->verts;
    if (pout->numverts < 3)
        return 0;

    return 1;
}

/////////////////////////////////////////////////////////////////////
// Clip a polygon to the frustum.
/////////////////////////////////////////////////////////////////////
static int ClipToFrustum(polygon_t* pin, polygon_t* pout) {
    int i, curpoly;
    polygon_t tpoly[2], *ppoly;

    curpoly = 0;
    ppoly = pin;

    for (i = 0; i < (NUM_FRUSTUM_PLANES - 1); i++) {
        if (!ClipToPlane(ppoly, &frustumplanes[i], &tpoly[curpoly])) {
            return 0;
        }
        ppoly = &tpoly[curpoly];
        curpoly ^= 1;
    }

    return ClipToPlane(ppoly, &frustumplanes[NUM_FRUSTUM_PLANES - 1], pout);
}

/////////////////////////////////////////////////////////////////////
// Add the polygon's edges to the global edge table.
/////////////////////////////////////////////////////////////////////
static void AddPolygonEdges(const plane_t* plane, polygon2D_t* screenpoly) {
    double distinv, deltax, deltay, slope;
    int i, nextvert, numverts, temp, topy, bottomy, height;
    edge_t* pedge;

    numverts = screenpoly->numverts;

    // Clamp the polygon's vertices just in case some very near
    // points have wandered out of range due to floating-point
    // imprecision
    for (i = 0; i < numverts; i++) {
        if (screenpoly->verts[i].x < -0.5)
            screenpoly->verts[i].x = -0.5;
        if (screenpoly->verts[i].x > (SCREEN_WIDTH - 0.5))
            screenpoly->verts[i].x = SCREEN_WIDTH - 0.5;
        if (screenpoly->verts[i].y < -0.5)
            screenpoly->verts[i].y = -0.5;
        if (screenpoly->verts[i].y > (SCREEN_HEIGHT - 0.5))
            screenpoly->verts[i].y = SCREEN_HEIGHT - 0.5;
    }

    // Add each edge in turn
    for (i = 0; i < numverts; i++) {
        nextvert = i + 1;
        if (nextvert >= numverts)
            nextvert = 0;

        topy = (int) ceil(screenpoly->verts[i].y);
        bottomy = (int) ceil(screenpoly->verts[nextvert].y);
        height = bottomy - topy;
        if (height == 0)
            continue; // doesn't cross any scan lines
        if (height < 0) {
            // Leading edge
            temp = topy;
            topy = bottomy;
            bottomy = temp;

            pavailedge->leading = 1;

            deltax = screenpoly->verts[i].x - screenpoly->verts[nextvert].x;
            deltay = screenpoly->verts[i].y - screenpoly->verts[nextvert].y;
            slope = deltax / deltay;

            // Edge coordinates are in 16.16 fixed point
            pavailedge->xstep = (int) (slope * (float) 0x10000);
            pavailedge->x =
                (int) ((screenpoly->verts[nextvert].x
                        + ((float) topy - screenpoly->verts[nextvert].y)
                              * slope)
                       * (float) 0x10000);
        } else {
            // Trailing edge
            pavailedge->leading = 0;

            deltax = screenpoly->verts[nextvert].x - screenpoly->verts[i].x;
            deltay = screenpoly->verts[nextvert].y - screenpoly->verts[i].y;
            slope = deltax / deltay;

            // Edge coordinates are in 16.16 fixed point
            pavailedge->xstep = (int) (slope * (float) 0x10000);
            pavailedge->x =
                (int) ((screenpoly->verts[i].x
                        + ((float) topy - screenpoly->verts[i].y) * slope)
                       * (float) 0x10000);
        }

        // Put the edge on the list to be added on top scan
        pedge = &newedges[topy];
        while (pedge->pnext->x < pavailedge->x)
            pedge = pedge->pnext;
        pavailedge->pnext = pedge->pnext;
        pedge->pnext = pavailedge;

        // Put the edge on the list to be removed after final scan
        pavailedge->pnextremove = removeedges[bottomy - 1];
        removeedges[bottomy - 1] = pavailedge;

        // Associate the edge with the surface we'll create for
        // this polygon
        pavailedge->psurf = pavailsurf;

        // Make sure we don't overflow the edge array
        if (pavailedge < &edges[MAX_EDGES])
            pavailedge++;
    }

    // Create the surface, so we'll know how to sort and draw from
    // the edges
    pavailsurf->state = 0;
    pavailsurf->color = currentcolor;

    // Set up the 1/z gradients from the polygon, calculating the
    // base value at screen coordinate 0,0 so we can use screen
    // coordinates directly when calculating 1/z from the gradients
    distinv = 1.0 / plane->distance;
    pavailsurf->zinvstepx =
        plane->normal.v[0] * distinv * maxscreenscaleinv * (fieldofview / 2.0);
    pavailsurf->zinvstepy =
        -plane->normal.v[1] * distinv * maxscreenscaleinv * (fieldofview / 2.0);
    pavailsurf->zinv00 = plane->normal.v[2] * distinv
                         - xcenter * pavailsurf->zinvstepx
                         - ycenter * pavailsurf->zinvstepy;

    // Make sure we don't overflow the surface array
    if (pavailsurf < &surfs[MAX_SURFS])
        pavailsurf++;
}

/////////////////////////////////////////////////////////////////////
// Scan all the edges in the global edge table into spans.
/////////////////////////////////////////////////////////////////////
static void ScanEdges(void) {
    int x, y;
    double fx, fy, zinv, zinv2;
    edge_t *pedge, *pedge2, *ptemp;
    span_t* pspan;
    surf_t *psurf, *psurf2;

    pspan = spans;

    // Set up the active edge list as initially empty, containing
    // only the sentinels (which are also the background fill). Most
    // of these fields could be set up just once at start-up
    edgehead.pnext = &edgetail;
    edgehead.pprev = NULL;
    edgehead.x = -0xFFFF; // left edge of screen
    edgehead.leading = 1;
    edgehead.psurf = &surfstack;

    edgetail.pnext = NULL; // mark edge of list
    edgetail.pprev = &edgehead;
    edgetail.x = SCREEN_WIDTH << 16; // right edge of screen
    edgetail.leading = 0;
    edgetail.psurf = &surfstack;

    // The background surface is the entire stack initially, and
    // is infinitely far away, so everything sorts in front of it.
    // This could be set just once at start-up
    surfstack.pnext = surfstack.pprev = &surfstack;
    surfstack.color = 0;
    surfstack.zinv00 = -999999.0;
    surfstack.zinvstepx = surfstack.zinvstepy = 0.0;

    for (y = 0; y < SCREEN_HEIGHT; y++) {
        fy = (double) y;

        // Sort in any edges that start on this scan
        pedge = newedges[y].pnext;
        pedge2 = &edgehead;
        while (pedge != &maxedge) {
            while (pedge->x > pedge2->pnext->x)
                pedge2 = pedge2->pnext;

            ptemp = pedge->pnext;
            pedge->pnext = pedge2->pnext;
            pedge->pprev = pedge2;
            pedge2->pnext->pprev = pedge;
            pedge2->pnext = pedge;

            pedge2 = pedge;
            pedge = ptemp;
        }

        // Scan out the active edges into spans

        // Start out with the left background edge already inserted,
        // and the surface stack containing only the background
        surfstack.state = 1;
        surfstack.visxstart = 0;

        for (pedge = edgehead.pnext; pedge; pedge = pedge->pnext) {
            psurf = pedge->psurf;

            if (pedge->leading) {
                // It's a leading edge. Figure out where it is
                // relative to the current surfaces and insert in
                // the surface stack; if it's on top, emit the span
                // for the current top.
                // First, make sure the edges don't cross
                if (++psurf->state == 1) {
                    fx = (double) pedge->x * (1.0 / (double) 0x10000);
                    // Calculate the surface's 1/z value at this pixel
                    zinv = psurf->zinv00 + psurf->zinvstepx * fx
                           + psurf->zinvstepy * fy;

                    // See if that makes it a new top surface
                    psurf2 = surfstack.pnext;
                    zinv2 = psurf2->zinv00 + psurf2->zinvstepx * fx
                            + psurf2->zinvstepy * fy;
                    if (zinv >= zinv2) {
                        // It's a new top surface
                        // emit the span for the current top
                        x = (pedge->x + 0xFFFF) >> 16;
                        pspan->count = x - psurf2->visxstart;
                        if (pspan->count > 0) {
                            pspan->y = y;
                            pspan->x = psurf2->visxstart;
                            pspan->color = psurf2->color;

                            // Make sure we don't overflow
                            // the span array
                            if (pspan < &spans[MAX_SPANS])
                                pspan++;
                        }

                        psurf->visxstart = x;

                        // Add the edge to the stack
                        psurf->pnext = psurf2;
                        psurf2->pprev = psurf;
                        surfstack.pnext = psurf;
                        psurf->pprev = &surfstack;
                    } else {
                        // Not a new top; sort into the surface stack.
                        // Guaranteed to terminate due to sentinel
                        // background surface
                        do {
                            psurf2 = psurf2->pnext;
                            zinv2 = psurf2->zinv00 + psurf2->zinvstepx * fx
                                    + psurf2->zinvstepy * fy;
                        } while (zinv < zinv2);

                        // Insert the surface into the stack
                        psurf->pnext = psurf2;
                        psurf->pprev = psurf2->pprev;
                        psurf2->pprev->pnext = psurf;
                        psurf2->pprev = psurf;
                    }
                }
            } else {
                // It's a trailing edge; if this was the top surface,
                // emit the span and remove it.
                // First, make sure the edges didn't cross
                if (--psurf->state == 0) {
                    if (surfstack.pnext == psurf) {
                        // It's on top, emit the span
                        x = ((pedge->x + 0xFFFF) >> 16);
                        pspan->count = x - psurf->visxstart;
                        if (pspan->count > 0) {
                            pspan->y = y;
                            pspan->x = psurf->visxstart;
                            pspan->color = psurf->color;

                            // Make sure we don't overflow
                            // the span array
                            if (pspan < &spans[MAX_SPANS])
                                pspan++;
                        }

                        psurf->pnext->visxstart = x;
                    }

                    // Remove the surface from the stack
                    psurf->pnext->pprev = psurf->pprev;
                    psurf->pprev->pnext = psurf->pnext;
                }
            }
        }

        // Remove edges that are done
        pedge = removeedges[y];
        while (pedge) {
            pedge->pprev->pnext = pedge->pnext;
            pedge->pnext->pprev = pedge->pprev;
            pedge = pedge->pnextremove;
        }

        // Step the remaining edges one scan line, and re-sort
        for (pedge = edgehead.pnext; pedge != &edgetail;) {
            ptemp = pedge->pnext;

            // Step the edge
            pedge->x += pedge->xstep;

            // Move the edge back to the proper sorted location,
            // if necessary
            while (pedge->x < pedge->pprev->x) {
                pedge2 = pedge->pprev;
                pedge2->pnext = pedge->pnext;
                pedge->pnext->pprev = pedge2;
                pedge2->pprev->pnext = pedge;
                pedge->pprev = pedge2->pprev;
                pedge->pnext = pedge2;
                pedge2->pprev = pedge;
            }

            pedge = ptemp;
        }
    }

    pspan->x = -1; // mark the end of the list
}

/////////////////////////////////////////////////////////////////////
// Draw all the spans that were scanned out.
/////////////////////////////////////////////////////////////////////
static void DrawSpans(void) {
    Uint8* buffer = VID_GetScreenBuffer();
    for (span_t* span = spans; span->x != -1; span++) {
        SDL_memset(buffer + (SCREEN_WIDTH * span->y) + span->x, span->color, span->count);
    }
}

/////////////////////////////////////////////////////////////////////
// Clear the lists of edges to add and remove on each scan line.
/////////////////////////////////////////////////////////////////////
static void ClearEdgeLists(void) {
    for (int i = 0; i < SCREEN_HEIGHT; i++) {
        newedges[i].pnext = &maxedge;
        removeedges[i] = NULL;
    }
}

/////////////////////////////////////////////////////////////////////
// Render the current state of the world to the screen.
/////////////////////////////////////////////////////////////////////
static void UpdateWorld(void) {
    polygon2D_t screenpoly;
    polygon_t *ppoly, tpoly0, tpoly1, tpoly2;
    convexobject_t* pobject;
    int i, j, k;
    plane_t plane;
    point_t tnormal;

    UpdateViewPos();
    SetUpFrustum();
    ClearEdgeLists();
    pavailsurf = surfs;
    pavailedge = edges;


    // Draw all visible faces in all objects
    pobject = objecthead.pnext;

    while (pobject != &objecthead) {
        ppoly = pobject->ppoly;

        for (i = 0; i < pobject->numpolys; i++) {
            // Move the polygon relative to the object center
            tpoly0.numverts = ppoly[i].numverts;
            for (j = 0; j < tpoly0.numverts; j++) {
                for (k = 0; k < 3; k++)
                    tpoly0.verts[j].v[k] =
                        ppoly[i].verts[j].v[k] + pobject->center.v[k];
            }

            if (PolyFacesViewer(&tpoly0, &ppoly[i].plane)) {
                if (ClipToFrustum(&tpoly0, &tpoly1)) {
                    currentcolor = ppoly[i].color;
                    TransformPolygon(&tpoly1, &tpoly2);
                    ProjectPolygon(&tpoly2, &screenpoly);

                    // Move the polygon's plane into viewspace
                    // First move it into worldspace (object relative)
                    tnormal = ppoly[i].plane.normal;
                    plane.distance = ppoly[i].plane.distance
                                     + DotProduct(&pobject->center, &tnormal);
                    // Now transform it into viewspace
                    // Determine the distance from the viewpont
                    plane.distance -= DotProduct(&currentpos, &tnormal);
                    // Rotate the normal into view orientation
                    plane.normal.v[0] = DotProduct(&tnormal, &vright);
                    plane.normal.v[1] = DotProduct(&tnormal, &vup);
                    plane.normal.v[2] = DotProduct(&tnormal, &vpn);

                    AddPolygonEdges(&plane, &screenpoly);
                }
            }
        }

        pobject = pobject->pnext;
    }

    ScanEdges();
    DrawSpans();
}


static void HandleKeyUp(const SDL_KeyboardEvent* e) {
    switch (e->keysym.scancode) {
        case SDL_SCANCODE_MINUS:
            fieldofview *= 0.9;
            xscreenscale = SCREEN_WIDTH / fieldofview;
            yscreenscale = SCREEN_HEIGHT / fieldofview;
            maxscale = SDL_max(xscreenscale, yscreenscale);
            break;
        case SDL_SCANCODE_EQUALS:
            fieldofview *= 1.1;
            xscreenscale = SCREEN_WIDTH / fieldofview;
            yscreenscale = SCREEN_HEIGHT / fieldofview;
            maxscale = SDL_max(xscreenscale, yscreenscale);
            break;
        case SDL_SCANCODE_F:
            speedscale *= 1.1;
            break;
        case SDL_SCANCODE_S:
            speedscale *= 0.9;
            break;
        default:
            break;
    }
}

static void HandleKeyDown(const SDL_KeyboardEvent* e) {
    switch (e->keysym.scancode) {
        case SDL_SCANCODE_DOWN:
            currentspeed -= MOVEMENT_SPEED * speedscale;
            currentspeed =
                SDL_max(currentspeed, -(MAX_MOVEMENT_SPEED * speedscale));
            break;
        case SDL_SCANCODE_UP:
            currentspeed += MOVEMENT_SPEED * speedscale;
            currentspeed =
                SDL_min(currentspeed, MAX_MOVEMENT_SPEED * speedscale);
            break;
        case SDL_SCANCODE_LEFT:
            yaw -= YAW_SPEED * speedscale;
            if (yaw < 0) {
                yaw += PI * 2;
            }
            break;

        case SDL_SCANCODE_RIGHT:
            yaw += YAW_SPEED * speedscale;
            if (yaw >= (PI * 2)) {
                yaw -= PI * 2;
            }
            break;
        case SDL_SCANCODE_N:
            roll += ROLL_SPEED * speedscale;
            if (roll >= (PI * 2)) {
                roll -= PI * 2;
            }
            break;
        case SDL_SCANCODE_M:
            roll -= ROLL_SPEED * speedscale;
            if (roll < 0) {
                roll += PI * 2;
            }
            break;
        case SDL_SCANCODE_A:
            pitch -= PITCH_SPEED * speedscale;
            if (pitch < 0) {
                pitch += PI * 2;
            }
            break;
        case SDL_SCANCODE_Z:
            pitch += PITCH_SPEED * speedscale;
            if (pitch >= (PI * 2)) {
                pitch -= PI * 2;
            }
            break;
        case SDL_SCANCODE_D:
            currentpos.v[1] += VMOVEMENT_SPEED;
            break;
        case SDL_SCANCODE_C:
            currentpos.v[1] -= VMOVEMENT_SPEED;
            break;
        default:
            break;
    }
}

void HandleEvents(void) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        switch (e.type) {
            case SDL_QUIT:
                quit = SDL_TRUE;
                break;
            case SDL_KEYDOWN:
                HandleKeyDown(&e.key);
                break;
            case SDL_KEYUP:
                HandleKeyUp(&e.key);
                break;
            default:
                break;
        }
    }
}

static SDL_bool InitInstance(void) {
    if (!VID_Init()) {
        printf("Could not initialize video! SDL_Error: %s\n", SDL_GetError());
        return SDL_FALSE;
    }

    // Set the initial location, direction, and speed
    roll = 0.0;
    pitch = 0.0;
    yaw = 0.0;
    currentspeed = 0.0;
    currentpos.v[0] = 0.0;
    currentpos.v[1] = 0.0;
    currentpos.v[2] = 0.0;
    fieldofview = 2.0;
    xscreenscale = SCREEN_WIDTH / fieldofview;
    yscreenscale = SCREEN_HEIGHT / fieldofview;
    maxscale = SDL_max(xscreenscale, yscreenscale);
    maxscreenscaleinv = 1.0 / maxscale;
    xcenter = SCREEN_WIDTH / 2.0 - 0.5;
    ycenter = SCREEN_HEIGHT / 2.0 - 0.5;

    numobjects = SDL_arraysize(objects);


    return SDL_TRUE; // We succeeded...
}

int main(int argc, char* argv[]) {
    if (!InitInstance()) {
        printf("Initialization failed!\n");
        return EXIT_FAILURE;
    }
    while (!quit) {
        HandleEvents();
        UpdateWorld();
        VID_UpdateScreen();
        SDL_Delay(16); // ~60 FPS
    }
    VID_Shutdown();
    return EXIT_SUCCESS;
}
