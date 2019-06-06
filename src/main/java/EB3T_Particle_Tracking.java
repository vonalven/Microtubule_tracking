import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

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

public class EB3T_Particle_Tracking implements PlugIn, DialogListener {

    private int nt;
    private ArrayList<Trajectory> trajectories = new ArrayList<Trajectory>();
    private Overlay overlay = new Overlay();
    private int nPasses = 1;                    // The number of passes (color channels * stack slices)
    private int pass;                           // Current pass

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

        // generate a dialog box to set the parameters for spots detection and trajectories creation
        GenericDialog dlg = new GenericDialog("Particle Tracking");
        dlg.addNumericField("Threshold", 0.4, 1);
        dlg.addNumericField("Radius spot [pixel]", 2.2, 1);
        dlg.addNumericField("Overlap intensity threshold (Holes threshold)", 2, 1);
        dlg.addNumericField("Velocity max [pixel/frame]", 10, 1);
        dlg.addCheckbox("Preview Detection", false);
        dlg.addDialogListener(this);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        double threshold = dlg.getNextNumber();
        double size = dlg.getNextNumber();
        double thresholdHole = dlg.getNextNumber();
        double velocity = dlg.getNextNumber();

        Thread thread = Thread.currentThread();     // needed to check for interrupted state
        long lastTime = System.currentTimeMillis();

        // 1) IDENTIFY ALL THE SPOTS FRAME-BY-FRAME
        // 2) BUILD TRAJECTORIES BY SPOTS LINKING
        for (int t = 1; t <= nt; t++) {
            imp.setPosition(t);
            ImageProcessor ip = imp.getProcessor();
            ImagePlus slice = new ImagePlus("", ip);
            ImagePlus dog = DoG(slice, size);
            ArrayList<Spot> spots = localMax(dog, t, threshold);
            for (Spot spot : spots) {
                spot.draw();
            }
            ArrayList<Spot> nonMatchedSpots = linkingNN(spots, velocity);
            for (Spot spot : nonMatchedSpots) {
                trajectories.add(new Trajectory(spot, trajectories.size() + 1));
            }

            long time = System.currentTimeMillis();
            if (time-lastTime > 100) {
                lastTime = time;
                if (thread.isInterrupted()) return;
                showProgress((double) Math.round(100-((double)t/(double)nt) * 100) / 100);
            }
        }


        ArrayList<Trajectory> newTrajectories = new ArrayList<>();

        // IF 2 TRAJECTORIES INTERSECATE, DETACH THE IDENTIFIED TRAJECTORIES IN THE CROSSING ZONE, THEY WILL BE BETTER
        // CONNECTED LATER. CROSSING IS IDENTIFIED BY A SHIFT-UP IN FLUORESCENCE INTENSITY
        // 1) IDENTIFY ALL THE SHIFT-UP IN INTENSITIES COMPARED TO THE AVERAGE INTENSITY OF THE TRAJECTORIES
        // 2) APPLY thresholdHole THRESHOLD TO IDENTIFY OVERLAP INTENSITIES PEAKS
        // 3) CUT THE TRAJECTORY IN THE IDENTIFIED CORSSING POINTS
        // 4) INSERT THE NEW TRAJECTORIES IN ALL THE TRAJECTORIES COLLECTION
        for (int idTraj = 0; idTraj < trajectories.size(); idTraj++) {
            Trajectory trajectory = trajectories.get(idTraj);
            double[] intensityDistribution = trajectory.getIntensityDistribution();
            for (int i = 0; i < intensityDistribution.length; i++) {
                if (intensityDistribution[i] > thresholdHole * trajectory.meanIntensity) {
                    intensityDistribution[i] = -1;
                }
            }

            int startPeakIndex = 0;
            int endPeakIndex = 0;
            int startSubTrajIndex = 0;
            boolean continuity = false;
            boolean hasPeak = false;

            for (int i = 0; i < intensityDistribution.length; i++) {
                if (intensityDistribution[i] == -1 && !continuity) {
                    startPeakIndex = i;
                    endPeakIndex = i;
                    continuity = true;
                    continue;
                }
                if (intensityDistribution[i] == -1 && continuity) {
                    endPeakIndex = i;
                    continue;
                }

                if (continuity && !(startPeakIndex == endPeakIndex)) {
                    if (startPeakIndex != 0) {
                        Trajectory tmpTraj = trajectory.getSubset(startSubTrajIndex, startPeakIndex); //TODO: verify
                        newTrajectories.add(tmpTraj);
                        hasPeak = true;
                    }

                    continuity = false;
                    startSubTrajIndex = endPeakIndex;
                }
            }

            if (hasPeak) {
                Trajectory tmpTraj = trajectory.getSubset(endPeakIndex, trajectory.size()); //TODO: verify
                newTrajectories.add(tmpTraj);

                trajectories.remove(trajectory);
            }
        }
        trajectories.addAll(newTrajectories);

        // draw the identified trajectories before linking. If the linking is performed this draw is replaced by the
        // linked trajectories
        for (Trajectory trajectory : trajectories) {
            trajectory.draw();
        }
        imp.setOverlay(overlay);

        // generate a dialog box to set the parameters for trajectories linking
        GenericDialog dlg2 = new GenericDialog("Match tracks");
        dlg2.addNumericField("Max gap [frames]", 10, 1);
        dlg2.addNumericField("Min sub-track length [frames]", 3, 1);
        dlg2.addNumericField("Forward cone half-opening [degrees]", 45, 1);
        dlg2.addCheckbox("Display statistics histograms", false);
        dlg2.addNumericField("Sedentary factor (higher removes sedentary points)", 0.25, 2);
        dlg2.showDialog();
        if (dlg2.wasCanceled()) return;
        double max_gap = dlg2.getNextNumber();
        double min_life = dlg2.getNextNumber();
        double fw_angle = dlg2.getNextNumber();
        boolean displayHisto = dlg2.getNextBoolean();

        double sedentaryFactor = dlg2.getNextNumber();

        ArrayList<Double> velocity_all = getVelocities(trajectories);

        double V_max = Percentile(velocity_all, 0.95);
        double V_med = Median(velocity_all);
        double bend_angle;
        double radius_max_fw;
        double delta_t_gap;
        double dist_test;
        Trajectory best_connected;
        boolean merge_event = true;
        int lastSize = trajectories.size();
        int counter = 0;

        // LINK TRAJECTORIES BY COST FUNCTION MINIMIZATION:
        // 1) IDENTIFY CANDIDATE LINKING TRAJECTORIES STARTING IN A FORWARD CONE OF HALF-APERTURE fw_angle AND A MAX
        //    RADIUS radius_max_fw = V_max * min(delta_t_gap, sqrt(max_gap))
        // 2) EXTRACT THE BEST CANDIDATE BY MINIMIZATION OF A COST FUNCTION
        // 3) LINK THE 2 TRAJECTORIES, REMOVE THE INDIVIDUAL BEST CANDIDATE FROM THE SET OF ALL TRAJECTORIES

        // repeat until nb. of trajectories (after linking) doesn't change
        while (merge_event) {

            long time = System.currentTimeMillis();
            if (time-lastTime > 100) {
                lastTime = time;
                if (thread.isInterrupted()) return;
                showProgress((double) Math.round(100-((double)counter/(double)trajectories.size()) * 100) / 100);
            }
            counter = 0; //just for the progress bar

            outerloop:
            for (int i = 0; i < trajectories.size(); i++) {
                counter++;
                Trajectory trajectory = trajectories.get(i);

                if (trajectory.isAlone) {
                    continue;
                }

                //TODO: da fare su ogni traiettoria PRIMA
                // build linear regression model of each Trajectory
                trajectory.buildRegression();

                // if the Trajectory is present on less than 2 frames, erase it (impossible to calculate bending and
                // lateral shift angles and an accurate linear regression model)
                if (trajectory.length_lin < 2) {
                    trajectories.remove(trajectory);
                    i--;
                    continue;
                }

                // If the Trajectory doesn't move a lot in the space is considered as sedentary.
                // Each sedentary trajectory will be removed.
                if (trajectory.isSedentary(velocity, sedentaryFactor)) {
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
                    best_connected = bestCandidate(trajectory, candidates_fw);

                    if (best_connected != null) {
                        trajectory.addAll(best_connected);
                        trajectories.remove(best_connected);
                        i--;
                    }
                } else {
                    trajectory.isAlone = true;
                }
            }

            if (lastSize == trajectories.size()) {
                merge_event = false;
            }
            lastSize = trajectories.size();
        }

        // remove trajectories shorter than min_life
        for (int i = 0; i < trajectories.size(); i++) {
            Trajectory trajectory = trajectories.get(i);
            if (trajectory.size() < min_life) {
                trajectories.remove(trajectory);
                i--;
            }
        }

        // display the trajectories
        ImagePlus imp2 = new ImagePlus();
        imp2 = imp.duplicate();
        imp2.show();
        //overlay = new Overlay();
        for (Trajectory trajectory : trajectories) {
            trajectory.draw();
        }
        imp2.setOverlay(overlay);

        // display the trajectories statistics in an imagej table
        ResultsTable table = new ResultsTable();
        for (Trajectory trajectory : trajectories) {
            table.incrementCounter();
            table.addValue("#", trajectory.num);
            table.addValue("length", trajectory.size());
            table.addValue("range", "(" + trajectory.start().t + " ... " + trajectory.last().t + ")");
            table.addValue("start position", "(" + trajectory.start().x + "," + trajectory.start().y + ")");
            table.addValue("last position", "(" + trajectory.last().x + "," + trajectory.last().y + ")");
        }
        table.show("Trajectories");

        // TODO
        // COMPUTE STATISTICS:
        ArrayList<Double> avg_velocities_all = getAvgVelocities(trajectories);
        ArrayList<Double> velocities_all = getVelocities(trajectories);
        ArrayList<Double> orientations_all = getOrientations(trajectories);
        ArrayList<Double> linear_lengths_all = getLinearLengths(trajectories);
        ArrayList<Double> curved_lengths_all = getCurvedLengths(trajectories);

        if (displayHisto) {
            Plot plot0 = new Plot("Average Velocity", "Average speed between spots for all trajectories [ pixels / frame ]", "Counts");
            double[] avg_V_arr = ArrayList_to_Array(avg_velocities_all);
            plot0.setColor(Color.red);
            plot0.setLineWidth(2);
            plot0.addHistogram(avg_V_arr);
            plot0.setLimitsToFit(true);
            plot0.show();

            Plot plot1 = new Plot("Velocity", "speed between spots [ pixels / frame ]", "Counts");
            double[] V_arr = ArrayList_to_Array(velocities_all);
            plot1.setColor(Color.red);
            plot1.setLineWidth(2);
            plot1.addHistogram(V_arr);
            plot1.setLimitsToFit(true);
            plot1.show();

            Plot plot2 = new Plot("Orientation", "Orientation angle between horizon and trajectories linear model [ degrees ]", "Counts");
            double[] or_arr = ArrayList_to_Array(orientations_all);
            plot2.setColor(Color.red);
            plot2.setLineWidth(2);
            plot2.addHistogram(or_arr);
            plot2.setLimitsToFit(true);
            plot2.show();

            Plot plot3 = new Plot("Linear Length", "Linear Length [ pixels ]", "Counts");
            double[] lin_len_arr = ArrayList_to_Array(linear_lengths_all);
            plot3.setColor(Color.red);
            plot3.setLineWidth(2);
            plot3.addHistogram(lin_len_arr);
            plot3.setLimitsToFit(true);
            plot3.show();

            Plot plot4 = new Plot("Curved Length", "Curved Length [ pixels ]", "Counts");
            double[] curv_len_arr = ArrayList_to_Array(curved_lengths_all);
            plot4.setColor(Color.red);
            plot4.setLineWidth(2);
            plot4.addHistogram(curv_len_arr);
            plot4.setLimitsToFit(true);
            plot4.setColor(Color.red);
            plot4.show();
        }

        GenericDialog dlg3 = new GenericDialog("Colorcale Plot");
        dlg3.addCheckbox("Avg velocity scale", false);
        dlg3.addCheckbox("Orientation scale", false);
        dlg3.addCheckbox("Linear length scale", false);
        dlg3.addCheckbox("Curved length scale", false);
        dlg3.showDialog();
        if (dlg3.wasCanceled()) return;
        boolean velocity_colorscale = dlg3.getNextBoolean();
        boolean orientation_colorscale = dlg3.getNextBoolean();
        boolean lin_length = dlg3.getNextBoolean();
        boolean curved_length = dlg3.getNextBoolean();

        if(velocity_colorscale) {
            ImagePlus imp_dupVelocity = new ImagePlus();
            imp_dupVelocity = imp.duplicate();
            imp_dupVelocity.setOverlay(overlay);
            drawWithColorGradient(trajectories, ArrayList_to_Array(avg_velocities_all), 100);
            imp_dupVelocity.show();
        }

        if(orientation_colorscale) {
            ImagePlus imp_dupOrientation = new ImagePlus();
            imp_dupOrientation = imp.duplicate();
            imp_dupOrientation.setOverlay(overlay);
            drawWithColorGradient(trajectories, ArrayList_to_Array(orientations_all), 100);
            imp_dupOrientation.show();
        }

        if(lin_length) {
            ImagePlus imp_dupLinLenght = new ImagePlus();
            imp_dupLinLenght = imp.duplicate();
            imp_dupLinLenght.setOverlay(overlay);
            drawWithColorGradient(trajectories, ArrayList_to_Array(linear_lengths_all), 100);
            imp_dupLinLenght.show();
        }

        if(curved_length) {
            ImagePlus imp_dupCurLenght = new ImagePlus();
            imp_dupCurLenght = imp.duplicate();
            imp_dupCurLenght.setOverlay(overlay);
            drawWithColorGradient(trajectories, ArrayList_to_Array(curved_lengths_all), 100);
            imp_dupCurLenght.show();
        }

    }

    private void drawWithColorGradient(ArrayList<Trajectory> trajectories, double[] statistics, int nb_colors) {
        int nb_col = nb_colors;
        if (nb_col < statistics.length) {
            nb_col = statistics.length;
        }

        Color[] colors = new Color[nb_col];
        double min_col = 0.0/360.0;
        double max_col = 120.0/360.0;
        double jump = (max_col - min_col) / (double) (nb_col);

        for (int i = 0; i < nb_col; i++) {
            colors[i] = Color.getHSBColor((float) (min_col + jump * i), 1.0f, 1.0f);
        }

        double[] ranges = range(getMin(statistics), getMax(statistics), nb_col);

        for (int j = 0; j<trajectories.size(); j++) {
            Trajectory traj = trajectories.get(j);
            loop_in:
            for (int i = 0; i < nb_col; i++) {
                if (i < nb_col - 2) {
                    if (statistics[j] >= ranges[i] && statistics[j] < ranges[i + 1]) {
                        traj.draw(colors[i]);
                        break loop_in;
                    }
                } else {
                    if (statistics[j] == ranges[i]) {
                        traj.draw(colors[nb_col - 1]);
                    }
                }
            }
        }
    }

    private double[] range(double min, double max, int nb_elements) {

        double[] arr = new double[nb_elements];
        double step = (max - min) / nb_elements;
        arr[0] = min;
        for (int i = 1; i < nb_elements - 1; i++) {
            arr[i] = arr[i - 1] + step;
        }
        arr[nb_elements - 1] = max;

        return arr;
    }

    public static double getMax(double[] inputArray) {
        double maxValue = inputArray[0];
        for (int i = 1; i < inputArray.length; i++) {
            if (inputArray[i] > maxValue) {
                maxValue = inputArray[i];
            }
        }
        return maxValue;
    }

    public static double getMin(double[] inputArray) {
        double minValue = inputArray[0];
        for (int i = 1; i < inputArray.length; i++) {
            if (inputArray[i] < minValue) {
                minValue = inputArray[i];
            }
        }
        return minValue;
    }

    // converts an ArrayList of Double to an array of double
    private double[] ArrayList_to_Array(ArrayList<Double> in) {
        double[] out = new double[in.size()];

        for (int i = 0; i < in.size(); i++) {
            out[i] = in.get(i);
        }

        return out;
    }


    private ArrayList<Double> getAvgVelocities(ArrayList<Trajectory> trajectories) {
        ArrayList<Double> avg_velocity_all = new ArrayList<>();

        double avg = 0;
        for (Trajectory trajectory : trajectories) {
            avg = 0;
            for (int i = 0; i < trajectory.size(); i++) {
                if (i < trajectory.size() - 1) {
                    avg += Math.abs(trajectory.get(i).distance(trajectory.get(i + 1)));
                }
            }
            avg_velocity_all.add(avg / (trajectory.size() - 1));
        }
        return avg_velocity_all;
    }

    // TODO
    // compute velocities [pixels/frame] bwtween all 2 successive points of all the identified trajectores. Store
    // them all in an unique array. Returns the array of velocities
    private ArrayList<Double> getVelocities(ArrayList<Trajectory> trajectories) {
        ArrayList<Double> velocity_all = new ArrayList<>();

        for (Trajectory trajectory : trajectories) {
            for (int i = 0; i < trajectory.size(); i++) {
                if (i < trajectory.size() - 1) {
                    velocity_all.add(Math.abs(trajectory.get(i).distance(trajectory.get(i + 1))));
                }
            }
        }
        return velocity_all;
    }

    // TODO
    // computes the orientations of all the trajectories in the input ArrayList. The orientations (angle) are computes
    // as atan(slope), where slope is given by the linear model attribute of each Trajectory
    private ArrayList<Double> getOrientations(ArrayList<Trajectory> trajectories) {
        ArrayList<Double> orientation_all = new ArrayList<>();

        for (Trajectory trajectory : trajectories) {
            orientation_all.add(Math.toDegrees(Math.atan(trajectory.regression.getSlope())));
        }
        return orientation_all;
    }

    // TODO
    // returns all the linear lengths of the trajectories in the input ArrayList. Those lengths are stored
    // as Trajectory class attributes
    private ArrayList<Double> getLinearLengths(ArrayList<Trajectory> trajectories) {
        ArrayList<Double> linear_lengths_all = new ArrayList<>();

        for (Trajectory trajectory : trajectories) {
            linear_lengths_all.add(trajectory.length_lin);
        }
        return linear_lengths_all;
    }

    // TODO
    // returns all the curved lengths of the trajectories in the input ArrayList. Those lengths are stored
    // as Trajectory class attributes
    private ArrayList<Double> getCurvedLengths(ArrayList<Trajectory> trajectories) {
        ArrayList<Double> curved_lengths_all = new ArrayList<>();

        for (Trajectory trajectory : trajectories) {
            curved_lengths_all.add(trajectory.length_curve);
        }
        return curved_lengths_all;
    }

    // return the element located at the percentile (<= 1) quantile of an input ArrayList
    private double Percentile(ArrayList<Double> ll, double percentile) {
        Collections.sort(ll);
        int index_percentile = (int) (Math.ceil(ll.size() * percentile));
        return ll.get(index_percentile);
    }

    // returns the median element of an input ArrayList
    private double Median(ArrayList<Double> ll) {
        int index_median = (int) (Math.ceil(ll.size() * 0.5));
        return ll.get(index_median);
    }

    // determines the best candidate linking to the reference Trajectory ref among a list of candidate Trajectories
    // (candidates_fw). The selection is based on the minimization of a cost function (cost_tmp_fw) considering the
    // lateral shift distance, the forward distance and the bending angle between the end and the begin of the old and
    // new candidate Trajectory
    private Trajectory bestCandidate(Trajectory ref,
                                     ArrayList<Trajectory> candidates_fw) {
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
                cost_tmp_fw = Math.abs(dist * Math.sin(bending_angle)) + dist * Math.abs(Math.cos(lateral_shift_angle)) + (1 - Math.cos(bending_angle));
                if (cost_tmp_fw < min_cost_fw) {
                    min_cost_fw = cost_tmp_fw;
                    best_cand_fw = cand_fw;
                }
            }
        }
        return best_cand_fw;
    }

    // computes the lateral shift angle between a reference Trajectory (ref) and a Spot (p). This angle corresponds to
    // the angle between the linear regression line of ref and the vector traced from the end of this regression to the
    // Spot p. It is equivalent to the half.angle of a cone with vertex lying on the "end of the linear regression" of
    // ref, with the outer aplitude-bonds defined by the Spot p.
    // isForward is an experimental parameter that was used to predict if the cone is oriented in the same growth
    // direction of ref or in the opposite direction (backward)
    private double shiftAngle(Trajectory ref, Spot p, boolean isForward) {
        double[] v1 = new double[2];
        double[] v2 = new double[2];
        v1[0] = p.x - ref.last().x;
        v1[1] = p.y - ref.regression.predict(ref.last().x);
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

    // override of shiftAngle. The angle is computed between the last spot of ref and the first spot of candidate
    // by calling the first shiftAngle method
    private double shiftAngle(Trajectory ref, Trajectory candidate, boolean isForward) {
        return shiftAngle(ref, candidate.start(), isForward);

    }

    // computes the bending angle between 2 Trajectory objects (ref = reference trajectory and candidate). This angle is
    // the angle between the 2 intersecating linear regression models specified as attributes in the Trajectory
    // objects
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

    // applies a Difference of Gaussians (DoG) filter on the input image in. The x and y sizes of the 2 gaussian sigmas
    // are respectively size and size * sqrt(2)
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

    // computes the local maximum based on 8-connected morphology of an input image (imp) using a grayscale threshold.
    // The found maxima are used to initialize an array of Spots.
    // t = frame of the imp processed image
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

    // links the detected Spots into trajectories by distance minimization. Link occurs if the minimal distance is
    // smaller than the velocity parameter and if the spot is in front of the growing closest trajectory. Otherwise the
    // spot is added in a nonMatchedSpots array
    public ArrayList<Spot> linkingNN(ArrayList<Spot> spots, double velocity) {
        ArrayList<Spot> nonMatchedSpots = new ArrayList<Spot>();
        for (Spot spot : spots) {
            Trajectory closest = null;
            double min = Double.MAX_VALUE;
            for (Trajectory trajectory : trajectories) {
                Spot p = trajectory.last();
                if (p.t == spot.t - 1) {
                    if (p.matched == null) {

                        double shiftAngle = 0;
                        if (trajectory.size() >= 2 && p.x != spot.x && p.y != spot.y) {
                            trajectory.buildRegression();
                            shiftAngle = shiftAngle(trajectory, spot, true) * unitPlateau(trajectory.size(), 0.5);
                        }


                        if (spot.distance(p) < min && shiftAngle <= 90) {
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

        private Color color;                    // color instance used by draw()
        private int num;                        // indentifier for the Trajectory number
        public SimpleRegression regression;     // linear regression model fitting the Spots of the Trajectory
        public double length_lin;               // "linar" length of the Trajectory
        public double length_curve;             // "curved" length of the Trajectory
        public boolean isAlone = false;         // True if the Trajectory is not linked to any other Trajectory
        public double[] intensityDistribution;  // intensity distribution at the Spots locations
        public double meanIntensity = 0;        // average intensity of intensityDistribution

        // constructor
        public Trajectory(Spot spot, int num) {
            this.num = num;
            if (spot.daughterOf != null) {
                add(spot.daughterOf.last());
                add(spot);
                color = spot.daughterOf.color;
            } else {
                add(spot);
                Random rand = new Random();
                // Java 'Color' class takes 3 floats, from 0 to 1.
                float r = rand.nextFloat();
                float g = rand.nextFloat();
                float b = rand.nextFloat();
                color = new Color(r, g, b, (float) 0.6);
            }
        }

        // constructor
        public Trajectory(List<Spot> subList) {
            for (int i = 0; i < subList.size(); i++) {
                this.add(subList.get(i));
            }
        }

        // returns the last Spot of the Trajectory
        public Spot last() {
            return get(this.size() - 1);
        }

        // returns the first Spot of the Trajectory
        public Spot start() {
            return get(0);
        }

        // draws all the points of the Trajectory and links them with lines
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

        // override of draw where a specific color can be specified
        public void draw(Color specific_color) {
            for (int i = 0; i < size() - 1; i++) {
                Spot a = get(i);
                Spot b = get(i + 1);
                Line line = new Line(a.x + 0.5, a.y + 0.5, b.x + 0.5, b.y + 0.5);
                line.setStrokeColor(specific_color);
                line.setStrokeWidth(1);
                overlay.add(line);
            }
        }

        // draws all the Spots of the Trajectory
        public void drawPoints() {
            for (int i = 0; i < size() - 1; i++) {
                get(i).draw();
            }
        }

        // draws the interpolated linear regression between the start and the end of the Trajectory
        public void drawRegression() {
            double ax = this.start().x;
            double ay = this.regression.predict(this.start().x);
            double bx = this.last().x;
            double by = this.regression.predict(this.last().x);
            Line line = new Line(ax, ay, bx, by);
            line.setStrokeColor(color);
            line.setStrokeWidth(1);
            overlay.add(line);
        }

        // computes a linear regression model fitting the Spots of the trajectory. Initializes the regression class
        // attribute
        // initializes the trajLenght_linear and the trajLength_curve class attributes
        public void buildRegression() {
            regression = new SimpleRegression();
            for (Spot spot : this) {
                regression.addData(spot.x, spot.y);
            }
            trajLenght_linear();
            trajLength_curve();
        }

        // computes the "linear" Trajectory length as a linear distance between the end and the start Spots evaluated in
        // the interpolated linear regression. Initializes the length_lin class attribute.
        public void trajLenght_linear() {
            this.length_lin = 0;
            double dx = this.last().x - this.start().x;
            double dy = this.regression.predict(this.last().x) - this.regression.predict(this.start().x);
            this.length_lin = Math.sqrt(dx * dx + dy * dy);
        }

        // computes the "curved" Trajectory length as sum of distances between all the Spots. Initializes the
        // length_curve class attribute.
        public void trajLength_curve() {
            this.length_curve = 0;
            for (int i = 0; i < this.size(); i++) {
                if (i < this.size() - 1) {
                    this.length_curve += Math.abs(this.get(i).distance(this.get(i + 1)));
                }
            }
        }

        // returns an array with the grayscale intensity of all the Spots in the Trajectory
        // initializes the intensityDistribution class attribute
        // computes the mean intensity of the intensityDistribution array
        public double[] getIntensityDistribution() {
            this.intensityDistribution = new double[this.size()];
            for (int i = 0; i < this.size(); i++) {
                this.intensityDistribution[i] = this.get(i).value;
            }
            this.meanIntensity = mean(this.intensityDistribution);
            return this.intensityDistribution;
        }

        // gets a subset of Spots from the current Trajectory
        public Trajectory getSubset(int startIndex, int endIndex) {
            return new Trajectory(this.subList(startIndex, endIndex));
        }

        public boolean isSedentary(double velocity, double factor) {
            double d = 0;
            for (int i = 1; i < this.size(); i++) {
                d += this.get(i).distance(this.get(0));
            }
            if ((d / this.size()) > (velocity * factor)) {
                return false;
            }
            return true;
        }
    }

    public class Spot {

        public int x;                   // x-coordinate of the Spot
        public int y;                   // y-coordinate of the Spot
        public int t;                   // frame at which the Spot is located
        public double value;            // grayscale intensity value of the input image at the Spot location
        public Spot matched;            // True if the Spot is linked to a trajectory
        public Trajectory daughterOf;   // trajectory linked to the Spot

        // constructor
        public Spot(int x, int y, int t, double value) {
            this.x = x;
            this.y = y;
            this.t = t;
            this.value = value;
        }

        // computes distance from current (.this) Spot and a second input spot
        public double distance(Spot spot) {
            double dx = x - spot.x;
            double dy = y - spot.y;
            return Math.sqrt(dx * dx + dy * dy);
        }

        // draw circles on Spot coordinates
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

    // displays a preview of the spots detected by DoG
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

    // exponential function used to scale values. f(x)=(1-e^(-a*x))
    // the plateu is reached at y = 1
    // a is the dilation factor
    private double unitPlateau(double x, double a) {
        return plateau(x, a, 1);
    }

    // exponential function used to scale values. f(x)=m*(1-e^(-a*x))
    // m is the plateau height (y_plateau = m)
    // a is a dilation factor
    private double plateau(double x, double a, double m) {
        return m * (1 - Math.exp(-a * x));
    }

    // computes average of an input array
    private double mean(double[] inputArray) {
        double sum = 0;

        for (int i = 0; i < inputArray.length; i++) {
            sum += inputArray[i];
        }

        int tot = inputArray.length != 0 ? inputArray.length : 1;
        return sum / tot;
    }


    /** This method is called by ImageJ to set the number of calls to run(ip)
     *  corresponding to 100% of the progress bar */
    public void setNPasses (int nPasses) {
        this.nPasses = nPasses;
        pass = 0;
    }

    private void showProgress(double percent) {
        percent = (double)(pass-1)/nPasses + percent/nPasses;
        IJ.showProgress(percent);
    }
}
