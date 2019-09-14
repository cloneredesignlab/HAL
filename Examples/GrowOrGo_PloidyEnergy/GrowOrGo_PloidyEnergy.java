package Examples.GrowOrGo_PloidyEnergy;

import Framework.GridsAndAgents.AgentGrid2D;
import Framework.GridsAndAgents.AgentSQ2Dunstackable;
import Framework.GridsAndAgents.PDEGrid2D;
import Framework.Gui.GridWindow;
import Framework.Tools.FileIO;
import Framework.Rand;
import Framework.Tools.MultiWellExperiment.MultiWellExperiment;
import Framework.Tools.PhylogenyTracker.Genome;
import Framework.Tools.PhylogenyTracker.GenomeFn;
import Framework.Util;

import java.io.FileReader;
import java.util.*;

import static Framework.Util.HeatMapRGB;
import static Framework.Util.WHITE;

class CopyNumber extends Genome {
    double[] cn_vals;
    static Rand rn = new Rand();
    int birthdate;

    public CopyNumber(Genome parent, double[] cn, int tick) {
        super(parent);
        this.cn_vals = cn.clone();
        this.birthdate = tick; //@TODO: birthdate tracking does not work as intended
    }

    CopyNumber[] Mutate(double mut_prob, int tick) {
        if (rn.Double() < mut_prob) {
            this.IncPop();
            this.IncPop();
            CopyNumber daughter1 = new CopyNumber(this, cn_vals, tick);
            CopyNumber daughter2 = new CopyNumber(this, cn_vals, tick);

            // 1st daughter
            int idx = rn.Int(cn_vals.length);
            if (daughter1.cn_vals[idx] != 0) { // Cell can't gain deletion back, nor can it loose additional copies
                double[] delta = new double[]{-1, 1}; //Loss, Gain
                daughter1.cn_vals[idx] += delta[rn.Int(2)];
            }

            //whatever is left of a perfect parental genome duplication goes to 2nd daughter
            for (int i = 0; i < daughter2.cn_vals.length; i++) {
                daughter2.cn_vals[i] = 2 * cn_vals[i] - daughter1.cn_vals[i];
            }
            return (new CopyNumber[]{daughter1, daughter2});
        } else {
            this.IncPop();
            return (new CopyNumber[]{this, this});
        }

    }

    double getPloidy() {
        double cnt = 0;
        for (double value : cn_vals) {
            cnt += value;
        }

//        System.out.println(cnt);
        return (cnt / (1.0 * cn_vals.length));
    }

}


//cells grow and mutate
class CellEx extends AgentSQ2Dunstackable<GrowOrGo_PloidyEnergy> {
    CopyNumber copy_number;
    double fitness;

    public CellEx() {
        super();
        double[] cn_vals = new double[22];
        Arrays.fill(cn_vals, G.PLOIDY);
        copy_number = new CopyNumber(null, cn_vals, 0);
        this.fitness = 1.0;
    }

    void Migrate(int nOpts) {
        G.energy.Add(Isq(), copy_number.getPloidy());

        int iNewLoc = G.hood[G.rn.Int(nOpts)];
        if (G.gradientMovement) {

            int x = G.ItoX(Isq());
            int y = G.ItoY(Isq());
            double[] grXY = new double[]{G.energy.GradientX(x, y), G.energy.GradientY(x, y)};

            //Calc. correlation oefficients
            Double[] corcoef = new Double[nOpts];
            for (int i = 0; i < nOpts; i++) {
                x = G.ItoX(G.hood[i]);
                y = G.ItoY(G.hood[i]);
                double[] grXY_ = new double[]{G.energy.GradientX(x, y), G.energy.GradientY(x, y)};
                corcoef[i] = correlation(grXY, grXY_);
                if (corcoef[i].equals(Double.NaN)) {
                    corcoef[i] = 0.0;
                }
            }
            double minCorr = Collections.min(Arrays.asList(corcoef));

            //Indices of maximum coefficients
            ArrayList<Integer> idx = new ArrayList<Integer>(nOpts);
            for (int i = 0; i < nOpts; i++) {
                if (corcoef[i] == minCorr) {
                    idx.add(G.hood[i]);
                }
            }
            iNewLoc = idx.get(G.rn.Int(idx.size()));
        }
        this.MoveSQ(iNewLoc);

        G.energy.Add(Isq(), copy_number.getPloidy() * -1);
    }

//    void Draw() {
//        G.vis.SetPix(Isq(), Util.HeatMapRGB(copy_number.getPloidy(), 1.7, 2.4));//sets a single pixel
//    }

    void Divide(int nOpts) {
        int iDaughter = G.hood[G.rn.Int(nOpts)];
        CopyNumber[] daughters = copy_number.Mutate(G.MUT_PROB, G.GetTick());//during division, there is a possibility of chromosome missegregation

        CellEx daughter = G.NewAgentSQ(iDaughter);//generate a daughter, the other is technically the original cell
        this.copy_number = daughters[0];
        daughter.copy_number = daughters[1];
        daughter.fitness = G.rn.Gaussian(this.fitness, 0.3); //Stochastic consequence of mutation on fitness
        this.fitness = G.rn.Gaussian(this.fitness, 0.3); //Stochastic consequence of mutation on fitness
        G.energy.Add(daughter.Isq(), daughter.copy_number.getPloidy() * -1);
    }

    void GrowOrGo(int nOpts, double availableEnergy, double ploidy) {
        if (availableEnergy >= Math.pow(ploidy, 3) && G.rn.Double() < G.DIV_PROB * fitness) {
            Divide(nOpts);
        } else if (G.rn.Double() < G.MIG_PROB) {
            Migrate(nOpts);
        }
    }

    void GoOrGrow(int nOpts, double availableEnergy, double ploidy) {
        if (G.rn.Double() < G.MIG_PROB) {
            Migrate(nOpts);
        } else if (availableEnergy >= Math.pow(ploidy, 3) && G.rn.Double() < G.DIV_PROB * fitness) {
            Divide(nOpts);
        }
    }

    //@TODO: move to utils
    public static double correlation(double[] xs, double[] ys) {
        double sx = 0.0;
        double sy = 0.0;
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;

        int n = xs.length;

        for (int i = 0; i < n; ++i) {
            double x = xs[i];
            double y = ys[i];

            sx += x;
            sy += y;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
        }

        // covariation
        double cov = sxy / n - sx * sy / n / n;
        // standard error of x
        double sigmax = Math.sqrt(sxx / n - sx * sx / n / n);
        // standard error of y
        double sigmay = Math.sqrt(syy / n - sy * sy / n / n);

        // correlation is just a normalized covariation
        return cov / sigmax / sigmay;
    }

}

public class GrowOrGo_PloidyEnergy extends AgentGrid2D<CellEx> {
    final static int BLACK = Util.RGB(0, 0, 0);
    double DIV_PROB = 0.22;
    double MUT_PROB = 0.02;
    double DIE_PROB = 0.15;
    double MIG_PROB = 0.3;
    double MUT_ADVANTAGE = 1.08;
    static double PLOIDY = 2.5;
    double DIFFUSION_RATE = 15;
    //    int[] hood = Util.GenHood2D(new int[]{1, 0, -1, 0, 0, 1, 0, -1}); //equivalent to int[]hood=Util.VonNeumannHood(false);
    int[] hood = Util.MooreHood(false);
    int[] hood_energy = Util.CircleHood(true, 4);
    boolean priority_grow;
    boolean gradientMovement;
    public PDEGrid2D energy;
    Rand rn = new Rand(22);
    FileIO outputFile = null;
    FileIO paramFile = null;

    public static void StepModel(GrowOrGo_PloidyEnergy model, int iWell) {
        model.StepCells();
    }

    public static int DrawModel(GrowOrGo_PloidyEnergy model, int x, int y) {
        CellEx c = model.GetAgent(x, y);
        if (c == null) {
            return BLACK; //Util.HeatMapJet(model.energy.Get(x, y));
        } else {
            return Util.HeatMapRGB(c.copy_number.getPloidy(), 2.1, 2.7);//sets a single pixel
        }
    }


    public GrowOrGo_PloidyEnergy(int x, int y, Map<String, String> argmap) {
        super(x, y, CellEx.class);
        this.priority_grow = Boolean.parseBoolean(argmap.get("priority_grow"));
        this.gradientMovement = Boolean.parseBoolean(argmap.get("gradientMovement"));
        this.DIE_PROB = Double.parseDouble(argmap.get("DIE_PROB"));
        this.MUT_PROB = Double.parseDouble(argmap.get("MUT_PROB"));
        this.MIG_PROB = Double.parseDouble(argmap.get("MIG_PROB"));
        this.MUT_ADVANTAGE = Double.parseDouble(argmap.get("MUT_ADVANTAGE"));
        PLOIDY = Double.parseDouble(argmap.get("PLOIDY"));
        this.DIFFUSION_RATE = Double.parseDouble(argmap.get("DIFFUSION_RATE"));
        energy = new PDEGrid2D(x, y);
        energy.SetAll(Double.parseDouble(argmap.get("energy")));
        energy.Update();
        this.priority_grow = Boolean.parseBoolean(argmap.get("priority_grow"));
        this.gradientMovement = Boolean.parseBoolean(argmap.get("gradientMovement"));
    }


    public GrowOrGo_PloidyEnergy(int x, int y, boolean priority_grow, boolean gradientMovement, double energyVal, String outputPrefix) {
        super(x, y, CellEx.class);
        this.priority_grow = priority_grow;
        this.gradientMovement = gradientMovement;
        if (outputPrefix != null) {
            String prfx = outputPrefix; // + "_priority_grow-" + priority_grow + "_energy" + energyVal;
            outputFile = new FileIO(prfx + ".out", "w");
            paramFile = new FileIO(prfx + ".param", "w");

            paramFile.Write("x\t" + x + "\n");
            paramFile.Write("y\t" + y + "\n");
            paramFile.Write("DIV_PROB\t" + DIV_PROB + "\n");
            paramFile.Write("MUT_PROB\t" + MUT_PROB + "\n");
            paramFile.Write("DIE_PROB\t" + DIE_PROB + "\n");
            paramFile.Write("MIG_PROB\t" + MIG_PROB + "\n");
            paramFile.Write("MUT_ADVANTAGE\t" + MUT_ADVANTAGE + "\n");
            paramFile.Write("PLOIDY\t" + PLOIDY + "\n");
            paramFile.Write("DIFFUSION_RATE\t" + DIFFUSION_RATE + "\n");
            paramFile.Write("energy\t" + energyVal + "\n");
            paramFile.Write("priority_grow\t" + priority_grow + "\n");
            paramFile.Write("gradientMovement\t" + gradientMovement);
            paramFile.Close();
        }
        energy = new PDEGrid2D(x, y);
        energy.SetAll(energyVal);
        energy.Update();
    }


    public void InitTumor(double radius) {
        //places tumor cells in a circle
        int[] circleHood = Util.CircleHood(true, radius);//generate circle neighborhood [x1,y1,x2,y2,...]
        int len = MapHood(circleHood, xDim / 2, yDim / 2);
        for (int i = 0; i < len; i++) {
            CellEx c = NewAgentSQ(circleHood[i]);
            energy.Add(c.Isq(), c.copy_number.getPloidy() * -1);
        }
        energy.Update();
    }

    public void StepCells() {
        //Cells
        for (CellEx c : this) {//iterate over all cells in the grid
            int nOpts = c.MapEmptyHood(hood);//finds von neumann neighborhood indices around cell.
            if (nOpts > 0) {
                double ploidy = c.copy_number.getPloidy();
//        int lockedEnergy = 2 * c.MapOccupiedHood(hood_energy);
//        double requiredEnergy = ploidy * 0.3 * c.MapHood(hood_energy);
                double availableEnergy = 0;
                for (int i = 0; i < c.MapEmptyHood(hood_energy); i++) {
                    availableEnergy += energy.Get(hood_energy[i]); //MAX_ENERGY - lockedEnergy - requiredEnergy;
                }
                availableEnergy = availableEnergy / hood_energy.length;

                if (rn.Double() < DIE_PROB / Math.pow(MUT_ADVANTAGE, ploidy)) { //application of ploidy advantage
                    c.Dispose();//removes cell from spatial grid and iteration
                    energy.Set(c.Isq(), ploidy); //Release cell's energy when it dies
                } else if (priority_grow) {
                    c.GrowOrGo(nOpts, availableEnergy, ploidy);
                } else {
                    c.GoOrGrow(nOpts, availableEnergy, ploidy);
                }
            }
        }
        ShuffleAgents(rn);//shuffles order of for loop iteration

        //Energy
        energy.DiffusionADI(DIFFUSION_RATE);
        energy.Update();
    }


    public static void mainShowSimulation(String paramFile) {

        ArrayList<String> lines = new FileIO(paramFile, "r").Read();
        Map<String, String> argmap = new HashMap<String, String>();

        for (String l : lines) {
            String[] kv = l.split("\t");
            argmap.put(kv[0], kv[1]);
        }

        int x = Integer.parseInt(argmap.get("x"));
        int y = Integer.parseInt(argmap.get("y"));
        GrowOrGo_PloidyEnergy[] model = new GrowOrGo_PloidyEnergy[]{new GrowOrGo_PloidyEnergy(x, y, argmap)};
        model[0].InitTumor(4);

        MultiWellExperiment<GrowOrGo_PloidyEnergy> expt = new MultiWellExperiment<GrowOrGo_PloidyEnergy>(1, 1, model, x, y,
                5, WHITE, GrowOrGo_PloidyEnergy::StepModel, GrowOrGo_PloidyEnergy::DrawModel, true);
        expt.RunGIF(50, paramFile + ".gif", 2, false);
//        expt.Run(50, false, 100);
    }


    public static void mainRunSimulations(int count, double eMin, double eMax, double by, String outDir) {
        int x = 70, y = 70, scaleFactor = 5;
        double[] energies = new double[(int) Math.round((eMax - eMin) / by)];
        for (int i = 0; i < energies.length; i++) {
            energies[i] = eMin + by * i;
        }

        int sim = 0;
        for (int n = 0; n < count; n++) {
            ////////////////////////
            ////Generate models////
            GrowOrGo_PloidyEnergy[] models = new GrowOrGo_PloidyEnergy[energies.length];
            double[] ploidies = new double[models.length];
            Arrays.fill(ploidies, 2);
            for (int i = 0; i < models.length; i++) {
                String lab = "yes";
                boolean priority_grow = true;
                if (i % 2 == 0) {
                    lab = "no";
                    priority_grow = false;
                }
                sim++; //Count this particular simulation
                models[i] = new GrowOrGo_PloidyEnergy(x, y, priority_grow, true, energies[i], outDir +
                        System.getProperty("file.separator") + "sim" + sim );
            }

            ////////////////////////
            ////Visualize models////
            for (int i = 0; i < models.length; i++) {
//      GrowOrGo_PloidyEnergy.PLOIDY = ploidies[i];
                models[i].InitTumor(4);
            }
            MultiWellExperiment<GrowOrGo_PloidyEnergy> expt = new MultiWellExperiment<GrowOrGo_PloidyEnergy>(6, 2, models, x, y,
                    scaleFactor, WHITE, GrowOrGo_PloidyEnergy::StepModel, GrowOrGo_PloidyEnergy::DrawModel, false);
            expt.Run(50, false, 100);


            ///////////////////////
            ////Analyze models////
            for (GrowOrGo_PloidyEnergy grid : models) {
                //Location and ploidy
                grid.outputFile.Write("Ploidy" + "\t" + "Dist2Center" + "\n");
                for (CellEx c : grid.AllAgents()) {
                    grid.outputFile.Write(c.copy_number.getPloidy() + "\t" + c.Dist(grid.xDim / 2, grid.yDim / 2) + "\n");
                }
                grid.outputFile.Close();

////      //Phylogeny @TODO: write output in Newick format
//                CopyNumber root = (CopyNumber) grid.AllAgents().iterator().next().copy_number.GetRoot();
//                String[] header = new String[22 + 3];
//                Arrays.fill(header, 0, 22, "chr");
//                header[22] = "BirthDate";
//                header[23] = "PopSize";
//                header[24] = "AliveSize";
//                grid.outputFile.Write(Util.ArrToString(header, "\t") + "\n");
//                root.Traverse(new GenomeFn() {
//                    @Override
//                    public void GenomeFn(Genome c) {
//                        grid.outputFile.Write(Util.ArrToString(((CopyNumber) c).cn_vals, "\t") + "\t" + c.GetId() + "\t" + c.GetNumGenomes() + "\t" + c.GetPop() + "\n");
//                    }
//                });
//                grid.outputFile.Close();
            }
        }
    }

    public static void main(String[] args) {
        Map<String, String> argmap = new HashMap<String, String>();
        for (String arg : args) {
            String[] kv = arg.split(":");
            argmap.put(kv[0], kv[1]);
        }

//                count:1 eMin:200 eMax:950 by:75 outDir:/Users/4470246/Downloads
        mainRunSimulations(Integer.parseInt(argmap.get("count")), Double.parseDouble(argmap.get("eMin")), Double.parseDouble(argmap.get("eMax")),
                Double.parseDouble(argmap.get("by")), argmap.get("outDir"));


//        //sim6_priority_grow-true_energy875.0.param
//        mainShowSimulation(args[0]);
    }

}