package Framework.GridsAndAgents;

import Framework.Interfaces.Coords2DToAction;
import Framework.Interfaces.Coords2DToBool;
import Framework.Util;

/**
 * Created by bravorr on 5/17/17.
 */
public abstract class GridBase2D extends GridBase{
    public final int xDim;
    public final int yDim;
    public final int length;
    public boolean wrapX;
    public boolean wrapY;
    public GridBase2D(int xDim, int yDim, boolean wrapX, boolean wrapY){
        this.xDim=xDim;
        this.yDim=yDim;
        this.wrapX=wrapX;
        this.wrapY=wrapY;
        this.length=xDim*yDim;
    }
    public int I(int x, int y){
        //gets typeGrid index from location
        return x*yDim+y;
    }

    /**
     * gets the index of the square at the specified coordinates with wrap around
     */
    public int WrapI(int x, int y){
        //wraps Coords to proper index
        if(In(x,y)) { return I(x,y);}
        return I(Util.ModWrap(x,xDim), Util.ModWrap(y,yDim));
    }

    /**
     * gets the xDim component of the square at the specified index
     */
    public int ItoX(int i){
        return i/yDim;
    }

    /**
     * gets the yDim component of the square at the specified index
     */
    public int ItoY(int i){
        return i%yDim;
    }

    /**
     * gets the index of the square that contains the specified coordinates
     */
    public int I(double x, double y){
        //gets typeGrid index from location
        return (int)Math.floor(x)*yDim+(int)Math.floor(y);
    }

    /**
     * returns whether the specified coordinates are inside the typeGrid bounds
     */
    public boolean In(int x, int y){
        return x >= 0 && x < xDim && y >= 0 && y < yDim;
    }

    /**
     * returns whether the specified coordinates are inside the typeGrid bounds
     */
    public boolean In(double x, double y){
        int xInt=(int)Math.floor(x);
        int yInt=(int)Math.floor(y);
        return In(xInt,yInt);
    }

//    /**
//     * gets the indices of the squares that lie within a given radius of a position
//     * argument array must be large enough to fit all indices in the maximum case, something like (rad*2)^2
//     */
//    public int SQsInRad(final int[] ret,final boolean wrapX,final boolean wrapY,final double centerX,final double centerY,final double rad){
//        int retCt=0;
//        for (int x = (int)Math.floor(centerX-rad); x <(int)Math.ceil(centerX+rad) ; x++) {
//            for (int y = (int)Math.floor(centerY-rad); y <(int)Math.ceil(centerY+rad) ; y++) {
//                int retX=x; int retY=y;
//                boolean inX=Util.InDim(xDim,retX);
//                boolean inY=Util.InDim(yDim,retY);
//                if((!wrapX&&!inX)||(!wrapY&&!inY)){
//                    continue;
//                }
//                if(wrapX&&!inX){
//                    retX=Util.ModWrap(retX,xDim);
//                }
//                if(wrapY&&!inY){
//                    retY=Util.ModWrap(retY,yDim);
//                }
//                ret[retCt]=I(retX,retY);
//                retCt++;
//            }
//        }
//        return retCt;
//    }
    /**
     * writes into ret the coordinates of SQs with the centerX and centerY coordinates added
     * returns the number of squares that were written into ret. any coordinates that fall out of bounds are not written
     * @param coords list of coordinates of the form [xDim,yDim,xDim,yDim,...]
     * @param ret list into which the displaced coordinates will be written
     * @param centerX xDim displacement of coordinates
     * @param centerY yDim displacement of coordinates
     * @param wrapX whether to wrap the coordinates that fall out of bounds in the X direction
     * @param wrapY whether to wrap the coordinates that fall out of bounds in the Y direction
     * @return the number of coordinates written into the ret array
     */
    public int HoodToIs(int[] coords, int[] ret, int centerX, int centerY, boolean wrapX, boolean wrapY){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2] + centerX;
            int y = coords[i * 2 + 1] + centerY;
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            ret[ptCt]= I(x,y);
            ptCt++;
        }
        return ptCt;
    }
    public int HoodToEvalIs(int[] coords,int []ret, int centerX, int centerY, Coords2DToBool Eval){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2] + centerX;
            int y = coords[i * 2 + 1] + centerY;
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            if(Eval.Eval(x,y)) {
                ret[ptCt]=I(x,y);
                ptCt++;
            }
        }
        return ptCt;
    }

    public int HoodToAction(int[] coords, int centerX, int centerY, Coords2DToAction Action){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2] + centerX;
            int y = coords[i * 2 + 1] + centerY;
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            Action.Action(x,y);
            ptCt++;
        }
        return ptCt;
    }
    public int CoordsToAction(int[] coords,boolean wrapX,boolean wrapY,Coords2DToAction Action){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2];
            int y = coords[i * 2 + 1];
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            Action.Action(x,y);
            ptCt++;
        }
        return ptCt;
    }
    public int HoodToIs(int[] coords, int[] ret, int centerX, int centerY){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2] + centerX;
            int y = coords[i * 2 + 1] + centerY;
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            ret[ptCt]= I(x,y);
            ptCt++;
        }
        return ptCt;
    }
    public int CoordsToIs(int[] coords,int[] ret,boolean wrapX,boolean wrapY){
            //moves coordinates to be around origin
            //if any of the coordinates are outside the bounds, they will not be added
            int ptCt=0;
            for(int i=0;i<coords.length/2;i++) {
                int x = coords[i * 2];
                int y = coords[i * 2 + 1];
                if (!Util.InDim(xDim, x)) {
                    if (wrapX) {
                        x = Util.ModWrap(x, xDim);
                    } else {
                        continue;
                    }
                }
                if (!Util.InDim(yDim, y)) {
                    if (wrapY) {
                        y = Util.ModWrap(y, yDim);
                    } else {
                        continue;
                    }
                }
                ret[ptCt]= I(x,y);
                ptCt++;
            }
            return ptCt;
    }
    public int CoordsToIs(int[] coords,int[] ret){
        //moves coordinates to be around origin
        //if any of the coordinates are outside the bounds, they will not be added
        int ptCt=0;
        for(int i=0;i<coords.length/2;i++) {
            int x = coords[i * 2];
            int y = coords[i * 2 + 1];
            if (!Util.InDim(xDim, x)) {
                if (wrapX) {
                    x = Util.ModWrap(x, xDim);
                } else {
                    continue;
                }
            }
            if (!Util.InDim(yDim, y)) {
                if (wrapY) {
                    y = Util.ModWrap(y, yDim);
                } else {
                    continue;
                }
            }
            ret[ptCt]= I(x,y);
            ptCt++;
        }
        return ptCt;
    }

    public double DistSquared(double x1, double y1, double x2, double y2, boolean wrapX, boolean wrapY){
        return Util.DistSquared(x1,y1,x2,y2, xDim, yDim, wrapX,wrapY);
    }
    public double DistSquared(double x1, double y1, double x2, double y2){
        return Util.DistSquared(x1,y1,x2,y2, xDim, yDim, wrapX,wrapY);
    }
    public int[] BoundaryIs(){
        int[] ret=new int[(xDim+yDim)*2];
        for (int x = 0; x < xDim; x++) {
            ret[x]=I(x,0);
            ret[x+xDim]=I(x,yDim-1);
        }
        for (int y = 0; y < yDim; y++) {
            ret[y+xDim*2]=I(0,y);
            ret[y+xDim*2+yDim]=I(xDim-1,y);
        }
        return ret;
    }
}
