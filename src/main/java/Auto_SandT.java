import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

public class Auto_SandT implements PlugIn {

    private int nt;
    private ArrayList<Trajectory> trajectories = new ArrayList<Trajectory>();
    private Overlay overlay = new Overlay();

    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("No open sequence of images.");
            return;
        }
        if (imp.getNFrames() < 1) {
            IJ.error("Requires a image with at least 2 frames.");
            return;
        }

        nt = imp.getStackSize();
        GenericDialog dlg = new GenericDialog("Auto S and T");
        dlg.addNumericField("Min Threshold", 0.5, 1);
        dlg.addNumericField("Max Threshold", 4.5, 1);
        dlg.addNumericField("Step Threshold", 0.100, 3);
        dlg.addNumericField("Min Size", 0.5, 1);
        dlg.addNumericField("Max Size", 4.5, 1);
        dlg.addNumericField("Step Size", 0.100, 3);
        dlg.addNumericField("Mean on n stacks", 2, 0);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        double minThreshold = dlg.getNextNumber();
        double maxThreshold = dlg.getNextNumber();
        double stepThreshold = dlg.getNextNumber();
        double minSize = dlg.getNextNumber();
        double maxSize = dlg.getNextNumber();
        double stepSize = dlg.getNextNumber();
        double tMean = dlg.getNextNumber();

        double[] thresholds = range(minThreshold,maxThreshold,stepThreshold);
        double sizes[] = range(minSize,maxSize,stepSize);

        int nThresholds = thresholds.length;
        int nSizes = sizes.length;

        int[] spotsEv = new int[nSizes];

        int threshold = 3;
        int velocity = 0;

        int randomStart = (int) (Math.random() * (nt-tMean) + 1); // numero da 0 a nt-tMean
        for(int s = 0; s <nSizes; s++) {
            double size = sizes[s];
            double nSpots = 0;
            for (int t = randomStart; t <= tMean; t++) {
                imp.setPosition(t);
                ImagePlus slice = new ImagePlus("", imp.getProcessor());
                ImagePlus dog = DoG(slice, size);
                ArrayList<Spot> spots = localMax(dog, t, threshold);
                nSpots += spots.size()/tMean;
            }
            spotsEv[s] = (int)Math.round(nSpots);
        }


        if(true){
            return;
        }



        for (Trajectory trajectory : trajectories) {
            trajectory.draw();
        }
        imp.setOverlay(overlay);

        ResultsTable table = new ResultsTable();
        for (Trajectory trajectory : trajectories) {
            table.incrementCounter();
            table.addValue("#", trajectory.num);
            table.addValue("length", trajectory.size());
            table.addValue("range", "(" + trajectory.start().t + " ... " + trajectory.last().t + ")");
            table.addValue("last position", "(" + trajectory.last().x + "," + trajectory.last().y + ")");
        }
        table.show("Trajectories");


    }

    public ImagePlus DoG(ImagePlus in, double size) {
        ImagePlus g1 = in.duplicate();
        ImagePlus g2 = in.duplicate();
        new ImageConverter(g1).convertToGray32();
        new ImageConverter(g2).convertToGray32();
        GaussianBlur3D.blur(g1, size, size, 0);
        GaussianBlur3D.blur(g2, size * Math.sqrt(2), size * Math.sqrt(2), 0);
        ImageProcessor ip1 = g1.getProcessor();
        ImageProcessor ip2 = g2.getProcessor();
        ip1.copyBits(ip2, 0, 0, Blitter.SUBTRACT);
        g1.setTitle("DoG of " + in.getTitle());
        return g1;
    }

    public ArrayList<Spot> localMax(ImagePlus imp, int t, double threshold) {
        int nx = imp.getWidth();
        int ny = imp.getHeight();
        ImageProcessor ip = imp.getProcessor();
        ArrayList<Spot> list = new ArrayList<Spot>();
        for (int x = 1; x < nx - 1; x++) {
            for (int y = 1; y < ny - 1; y++) {
                double v = ip.getPixelValue(x, y);
                if (v > threshold) {
                    if (v > ip.getPixelValue(x - 1, y + 1)) {
                        if (v > ip.getPixelValue(x - 1, y)) {
                            if (v > ip.getPixelValue(x - 1, y - 1)) {
                                if (v > ip.getPixelValue(x + 1, y + 1)) {
                                    if (v > ip.getPixelValue(x + 1, y)) {
                                        if (v > ip.getPixelValue(x + 1, y - 1)) {
                                            if (v > ip.getPixelValue(x, y - 1)) {
                                                if (v > ip.getPixelValue(x, y + 1)) {
                                                    list.add(new Spot(x, y, t, v));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return list;
    }

    public ArrayList<Spot> linkingNN(ArrayList<Spot> spots, double velocity) {
        ArrayList<Spot> nonMatchedSpots = new ArrayList<Spot>();
        for (Spot spot : spots) {
            Trajectory closest = null;
            double min = Double.MAX_VALUE;
            for (Trajectory trajectory : trajectories) {
                Spot p = trajectory.last();
                if (p.t == spot.t - 1)
                    if (p.matched == null)
                        if (spot.distance(p) < min) {
                            closest = trajectory;
                            min = spot.distance(p);
                        }
            }
            if (min < velocity) {
                closest.last().matched = spot;
                closest.add(spot);
            }
            else {
                nonMatchedSpots.add(spot);
            }
        }
        return nonMatchedSpots;
    }

    public class Trajectory extends ArrayList<Spot> {

        private Color color;
        private int num;

        public Trajectory(Spot spot, int num) {
            this.num = num;
            if (spot.daughterOf == null) {
                add(spot);
                color = Color.getHSBColor((float) Math.random(), 1f, 1f);
                color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 150);
            }
            else {
                add(spot.daughterOf.last());
                add(spot);
                color = spot.daughterOf.color;
            }
        }

        public Spot last() {
            return get(size() - 1);
        }

        public Spot start() {
            return get(0);
        }

        public void draw() {
            //overlay.add(new TextRoi(last().x, last().y, ""+num));
            for (int i = 0; i < size() - 1; i++) {
                Spot a = get(i);
                Spot b = get(i + 1);
                Line line = new Line(a.x + 0.5, a.y + 0.5, b.x + 0.5, b.y + 0.5);
                line.setStrokeColor(color);
                line.setStrokeWidth(1);
                overlay.add(line);
            }
        }
    }

    public class Spot {

        public int x;
        public int y;
        public int t;
        public double value;
        public Spot matched;
        public Trajectory daughterOf;

        public Spot(int x, int y, int t, double value) {
            this.x = x;
            this.y = y;
            this.t = t;
            this.value = value;
        }

        public double distance(Spot spot) {
            double dx = x - spot.x;
            double dy = y - spot.y;
            return Math.sqrt(dx * dx + dy * dy);
        }

        public void draw() {
            double xp = x + 0.5;
            double yp = y + 0.5;
            int radius = 5;
            Roi roi = new OvalRoi(xp - radius, yp - radius, 2 * radius, 2 * radius);
            roi.setPosition(t);
            roi.setStrokeColor(Color.WHITE);
            roi.setStrokeWidth(1);
            overlay.add(roi);
        }
    }

    private double[] range(double min, double max, double step){
        int nEls = (int)(Math.ceil((max-min)/step))+1;
        double arr[] = new double[nEls];
        arr[0] = min;
        for(int i = 1; i < nEls-1; i++){
            arr[i] = Math.round((arr[i-1]+step)*1000.0)/1000.0;
            IJ.log(String.valueOf(arr[i]));
        }
        arr[nEls-1] = max;

        return arr;
    }
}
