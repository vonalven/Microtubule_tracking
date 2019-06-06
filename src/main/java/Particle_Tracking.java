import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import org.apache.commons.math3.stat.regression.SimpleRegression;

public class Particle_Tracking implements PlugIn, DialogListener {

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
        imp.setOverlay(overlay);

        nt = imp.getStackSize();
        GenericDialog dlg = new GenericDialog("Particle Tracking");
        dlg.addNumericField("Threshold", 10, 1);
        dlg.addNumericField("Radius spot [pixel]", 3, 1);
        dlg.addNumericField("Velocity max [pixel/frame]", 10, 1);
        dlg.addCheckbox("Preview Detection", false);
        dlg.addDialogListener(this);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        double threshold = dlg.getNextNumber();
        double size = dlg.getNextNumber();
        double velocity = dlg.getNextNumber();

        for (int t = 1; t <= nt; t++) {
            imp.setPosition(t);
            ImagePlus slice = new ImagePlus("", imp.getProcessor());
            ImagePlus dog = DoG(slice, size);
            ArrayList<Spot> spots = localMax(dog, t, threshold);
            for (Spot spot : spots) {
                spot.draw();
            }
            ArrayList<Spot> nonMatchedSpots = linkingNN(spots, velocity);
            for (Spot spot : nonMatchedSpots) {
                trajectories.add(new Trajectory(spot, trajectories.size() + 1));
            }
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

        GenericDialog dlg2 = new GenericDialog("Match tracks");
        dlg2.addNumericField("Max gap [frames]", 10, 1);
        dlg2.addNumericField("Min sub-track length [frames]", 3, 1);
        dlg2.addNumericField("Forward cone half-opening [degrees]", 45, 1);
        dlg2.addNumericField("Backward cone half-opening [degrees]", 10, 1);
        dlg2.addNumericField("Max transveral shift angle [degrees]", 20, 1);
        dlg2.addNumericField("Backward velocity shrinkage coeff", 1.5, 2);
        dlg2.addDialogListener(this);
        dlg2.showDialog();
        if (dlg.wasCanceled()) return;
        double max_gap = dlg2.getNextNumber();
        double min_life = dlg2.getNextNumber();
        double fw_angle = dlg2.getNextNumber();
        double bkw_angle = dlg2.getNextNumber();
        double transversal_angle = dlg2.getNextNumber();
        double V_shrinkage = dlg2.getNextNumber();
        ArrayList<Double> velocity_all = new ArrayList<>();

        for (Trajectory trajectory : trajectories) {
            for (int i = 0; i < trajectory.size(); i++) {
                if (i < trajectory.size() - 1) {
                    velocity_all.add(Math.abs(trajectory.get(i).distance(trajectory.get(i + 1))));
                }
            }
        }

        double V_max = Percentile(velocity_all, 0.95);
        double V_med = Median(velocity_all);
        double shift_angle_fw;
        double bend_angle;
        double radius_max_fw;
        double delta_t_gap;
        double dist_test;
        Trajectory best_connected;
        boolean merge_event = true;

        // repeat until nb. of trajectories (after linking) doesn't change
        while (merge_event) {

            merge_event = false;

            outerloop:
            for (int i = 0; i < trajectories.size(); i++) {
                Trajectory trajectory = trajectories.get(i);

                if (trajectory.isAlone) {
                    continue;
                }

                //TODO: da fare su ogni traiettoria PRIMA
                trajectory.buildRegression();

                if(trajectory.length_lin < 2){
                    trajectories.remove(trajectory);
                    i--;
                    continue;
                }

                if(trajectory.length_curve > trajectory.length_lin*1.5){
                    trajectories.remove(trajectory);
                    i--;
                    continue;
                }

                ArrayList<Trajectory> candidates_fw = new ArrayList<>();

                for (Trajectory traj : trajectories) {
                    if (traj != trajectory) {

                        traj.buildRegression();
                        delta_t_gap = Math.abs(trajectory.last().t - traj.start().t);

                        if (delta_t_gap < max_gap) {

                            bend_angle = bendAngle(trajectory, traj);
                            double bend_corrected = bend_angle * Math.tanh(traj.length_lin / max_gap);

                            radius_max_fw = V_max * Math.min(delta_t_gap, Math.sqrt(max_gap));

                            dist_test = Math.abs(trajectory.last().distance((traj.start())));

                            if (dist_test <= radius_max_fw && bend_corrected <= fw_angle) {
                                candidates_fw.add(traj);
                            }
                        }
                    }
                }

                if (!candidates_fw.isEmpty()) {
                    best_connected = bestCandidate(trajectory, candidates_fw, fw_angle, transversal_angle);

                    if(best_connected!=null) {
                        trajectory.addAll(best_connected);
                        trajectories.remove(best_connected);
                        merge_event = true;
                    }
                    break outerloop;
                } else {
                    trajectory.isAlone = true;
                }
            }
        }

        for (int i = 0; i < trajectories.size(); i++) {
            Trajectory trajectory = trajectories.get(i);
            if(trajectory.size() < min_life){
                trajectories.remove(trajectory);
                i--;
            }
        }


        ImagePlus imp2 = new ImagePlus();
        imp2 = imp.duplicate();
        imp2.show();
        //overlay = new Overlay();
        for (Trajectory trajectory : trajectories) {
            trajectory.draw();
        }
        imp2.setOverlay(overlay);
    }

    private double Percentile(ArrayList<Double> ll, double percentile) {
        Collections.sort(ll);
        int index_percentile = (int) (Math.ceil(ll.size() * percentile));
        return ll.get(index_percentile);
    }

    private double Median(ArrayList<Double> ll) {
        int index_median = (int) (Math.ceil(ll.size() * 0.5));
        return ll.get(index_median);
    }

    // TODO: last 2 parameters never used
    private Trajectory bestCandidate(Trajectory ref,
                                     ArrayList<Trajectory> candidates_fw,
                                     double max_angle_fw, double max_angle_lateral) {
        double bending_angle;
        double lateral_shift_angle;
        double cost_tmp_fw;
        double min_cost_fw = Double.MAX_VALUE;
        double dist;
        Trajectory best_cand_fw = null;

        if (!candidates_fw.isEmpty()) {
            for (Trajectory cand_fw : candidates_fw) {
                lateral_shift_angle = shiftAngle(ref, cand_fw, true);
                bending_angle = bendAngle(ref, cand_fw);

                dist = Math.abs(ref.last().distance(cand_fw.start()));
                //cost_tmp_fw = Math.abs(Math.cos(bending_angle)) - Math.abs(Math.cos(lateral_shift_angle));
                // cost_tmp_fw = Math.abs(Math.sin(bending_angle)) + Math.abs(Math.cos(lateral_shift_angle));
                cost_tmp_fw = Math.abs(dist*Math.sin(bending_angle)) + dist*Math.abs(Math.cos(lateral_shift_angle)) + (1-Math.cos(bending_angle));
                if (cost_tmp_fw < min_cost_fw) {
                    min_cost_fw = cost_tmp_fw;
                    best_cand_fw = cand_fw;
                }
            }
        }
        return best_cand_fw;
    }

    private double shiftAngle(Trajectory ref, Trajectory candidate, boolean isForward) {
        double[] v1 = new double[2];
        double[] v2 = new double[2];
        v1[0] = candidate.start().x - ref.last().x;
        v1[1] = candidate.start().y - ref.regression.predict(ref.last().x);
        if (isForward) {
            v2[0] = ref.last().x - ref.start().x;
            v2[1] = ref.regression.predict(ref.last().x) - ref.regression.predict(ref.start().x);
        }

        double norm1 = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
        double norm2 = Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
        v1[0] /= norm1;
        v1[1] /= norm1;
        v2[0] /= norm2;
        v2[1] /= norm2;

        return Math.abs(Math.toDegrees(Math.acos(v1[0] * v2[0] + v1[1] * v2[1])));

    }

    private double bendAngle(Trajectory ref, Trajectory candidate) {

        double[] v1 = new double[2];
        double[] v2 = new double[2];
        v1[0] = ref.start().x - ref.last().x;
        v1[1] = ref.regression.predict(ref.start().x) - ref.regression.predict(ref.last().x);

        v2[0] = candidate.start().x - candidate.last().x;
        v2[1] = candidate.regression.predict(ref.start().x) - candidate.regression.predict(ref.last().x);

        double norm1 = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
        double norm2 = Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
        v1[0] /= norm1;
        v1[1] /= norm1;
        v2[0] /= norm2;
        v2[1] /= norm2;

        return Math.abs(Math.toDegrees(Math.acos(v1[0] * v2[0] + v1[1] * v2[1])));

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
                if (p.t == spot.t - 1) {
                    if (p.matched == null) {
                        if (spot.distance(p) < min) {
                            closest = trajectory;
                            min = spot.distance(p);
                        }
                    }
                }
            }
            if (!(min < velocity)) {
                nonMatchedSpots.add(spot);
            } else {
                closest.last().matched = spot;
                closest.add(spot);
            }
        }
        return nonMatchedSpots;
    }

    public class Trajectory extends ArrayList<Spot> {

        private Color color;
        private int num;
        public SimpleRegression regression;
        public double length_lin;
        public double length_curve;
        public boolean isAlone = false;

        public Trajectory(Spot spot, int num) {
            this.num = num;
            if (spot.daughterOf != null) {
                add(spot.daughterOf.last());
                add(spot);
                color = spot.daughterOf.color;
            } else {
                add(spot);
                color = Color.getHSBColor((float) Math.random(), 1f, 1f);
                color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 150);
            }
        }

        public Spot last() {
            return get(this.size() - 1);
        }

        public Spot start() {
            return get(0);
        }

        public void draw() {
            for (int i = 0; i < size() - 1; i++) {
                Spot a = get(i);
                Spot b = get(i + 1);
                Line line = new Line(a.x + 0.5, a.y + 0.5, b.x + 0.5, b.y + 0.5);
                line.setStrokeColor(color);
                line.setStrokeWidth(1);
                overlay.add(line);
            }
        }

        public void buildRegression() {
            regression = new SimpleRegression();
            for (Spot spot : this) {
                regression.addData(spot.x, spot.y);
            }
            trajLenght_linear();
            trajLength_curve();
        }

        public void trajLenght_linear() {
            this.length_lin = 0;
            double dx = this.last().x - this.start().x;
            double dy = this.regression.predict(this.last().x) - this.regression.predict(this.start().x);
            this.length_lin = Math.sqrt(dx * dx + dy * dy);
        }

        public void trajLength_curve(){
            this.length_curve = 0;
            for (int i = 0; i<this.size(); i++) {
                if (i<this.size()-1){
                    this.length_curve += Math.abs(this.get(i).distance(this.get(i+1)));
                }
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

    @Override
    public boolean dialogItemChanged(GenericDialog dlg, AWTEvent event) {
        double threshold = dlg.getNextNumber();
        double size = dlg.getNextNumber();
        ImagePlus imp = IJ.getImage();
        overlay = new Overlay();
        if (dlg.getNextBoolean()) {
            ImagePlus slice = new ImagePlus("", imp.getProcessor());
            int t = imp.getFrame();
            ImagePlus dog = DoG(slice, size);
            ArrayList<Spot> spots = localMax(dog, t, threshold);
            for (Spot spot : spots)
                spot.draw();
        }
        imp.setOverlay(overlay);
        return !dlg.invalidNumber();
    }
}
