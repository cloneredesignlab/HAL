package Testing.OldTests;

import Framework.GridsAndAgents.PDEGrid1D;
import Framework.GridsAndAgents.PDEGrid2D;
import Framework.GridsAndAgents.PDEGrid3D;
import Framework.Gui.GridWindow;
import Framework.Util;

import static Framework.Util.HeatMapRGB;

public class AdvectionTest1D {
    public static void main(String[] args) {
        GridWindow win=new GridWindow(100,10,10);
        PDEGrid3D grid=new PDEGrid3D(100,10,10);
        while(true){
            win.TickPause(0);
            grid.Advection(0.1,0.1,0.1);
            grid.Set(0,1);
            grid.Update();
            win.DrawPDEGridXY(grid, Util::HeatMapRGB);
        }
    }
}
