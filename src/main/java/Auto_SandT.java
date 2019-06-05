import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.GaussianBlur3D;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

public class Auto_SandT implements PlugIn {

    private int nt;
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

        double veryMin = 0.5;
        double minIntensity = Math.max(imp.getProcessor().getMin(),veryMin);
        double autoIntensity = imp.getProcessor().getAutoThreshold();

        if(minIntensity>autoIntensity){
            double tmp = autoIntensity;
            autoIntensity = minIntensity;
            minIntensity = tmp;
        }

        IJ.log(String.valueOf(autoIntensity));

        double threshold = autoIntensity;//(autoIntensity+minIntensity)/2;

        //TODO: ask to the user with a dialog (see old version)
        double minSize = 0.2;
        double maxSize = 4.5;
        double stepSize = 0.1;
        double tMean = 4;

        double sizes[] = range(minSize,maxSize,stepSize);

        int nSizes = sizes.length;


        int velocity = 0;

        int randomStart = (int) (Math.random() * (nt-tMean) + 1); // numero da 0 a nt-tMean

        double[] spotsEv = new double[nSizes];
        for(int s = 0; s <nSizes; s++) {
            double size = sizes[s];
            double nSpots = 0;
            for (int t = randomStart; t < randomStart+tMean; t++) {
                imp.setPosition(t);
                ImagePlus slice = new ImagePlus("", imp.getProcessor());
                ImagePlus dog = DoG(slice, size);
                ArrayList<Spot> spots = localMax(dog, t, threshold);
                nSpots += spots.size()/tMean;
            }
            spotsEv[s] = Math.round(nSpots);
        }

        Plot plot = new Plot("prova th: " + threshold, "size", "counts");
        plot.addPoints(sizes, spotsEv, Plot.LINE);
        plot.setColor(Color.green);

        double[] deriv    = derivative(spotsEv,5,true);
        //normalize:
        deriv = normalize(deriv);
        double[] newSizes = new double[deriv.length];
        for(int fac = 0; fac<deriv.length;fac++){
            newSizes[fac] = sizes[fac+2];
        }

        plot.addPoints(newSizes,deriv,Plot.BOX);
        plot.setLimitsToFit(true);
        plot.show();


        int[] bornes = stableMoovingMean(deriv);
        double bestSize = (sizes[bornes[0]]-sizes[bornes[1]])/2;
        IJ.log("best size: "+bestSize);


        if(true){
            return;
        }

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

    public class Spot {

        public int x;
        public int y;
        public int t;
        public double value;
        public Spot matched;

        public Spot(int x, int y, int t, double value) {
            this.x = x;
            this.y = y;
            this.t = t;
            this.value = value;
        }
    }

    private double[] derivative2(double[] original,int precision,boolean abs){
        double[] d1 = derivative(original,precision,false);
        return derivative(d1,precision,abs);
    }

    private double[] derivative(double[] original){
        return derivative(original,2,false);
    }

    private double[] derivativeAbsolute(double[] original){
        return derivative(original,2,true);
    }

    private double[] derivativeAbsolute(double[] original,int precision){
        return derivative(original,precision,true);
    }

    private double[] derivative(double[] original,int precision, boolean abs){
        int l = original.length - precision*2;
        double[] d = new double[l];
        for(int x=precision;x<l-precision;x++){
            d[x-precision] = (original[x+precision] - original[x-precision]) / (2*precision);
            if(abs){
                d[x-precision] = Math.abs(d[x-precision]);
            }
        }
        return d;
    }

    private double[] normalize(double[] original){
        double max = getMax(original);
        double[] outArray = new double[original.length];
        for (int i = 0; i < original.length; i++) {
            outArray[i] = original[i]/max;
        }
        return outArray;
    }

    private double[] range(double min, double max, double step){
        int nEls = (int)(Math.ceil((max-min)/step))+1;
        double arr[] = new double[nEls];
        arr[0] = min;
        for(int i = 1; i < nEls-1; i++){
            arr[i] = Math.round((arr[i-1]+step)*1000.0)/1000.0;
        }
        arr[nEls-1] = max;

        return arr;
    }

    // Method for getting the maximum value
    public static double getMax(double[] inputArray){
        double maxValue = inputArray[0];
        for(int i=1;i < inputArray.length;i++){
            if(inputArray[i] > maxValue){
                maxValue = inputArray[i];
            }
        }
        return maxValue;
    }

    private int[] stableMoovingMean(double[] inputArray){
        return stableMoovingMean(inputArray,0.2,10,20);
    }

    private int[] stableMoovingMean(double[] inputArray, double threshold, int window, int limit){
        int[] minAndMax = new int[2];
        double mean = 0;
        double lastMean = 0;
        boolean checkpoint = false;
        for (int i = 0; i < (inputArray.length-window); i++) {
            for (int j = 0; j < window; j++) {
                mean += inputArray[i+j];
            }
            mean = mean/window;
            if(!checkpoint && Math.abs(mean-lastMean)<threshold){
                minAndMax[0] = i;
                checkpoint = true;
                continue;
            }
            if(checkpoint && Math.abs(mean-lastMean)>threshold){
                minAndMax[1] = i;
                break;
            }
            if(i>limit){
                minAndMax[1] = i;
                break;
            }
        }
        return minAndMax;
    }

    // Method for getting the minimum value
    public static double getMin(double[] inputArray){
        double minValue = inputArray[0];
        for(int i=1;i<inputArray.length;i++){
            if(inputArray[i] < minValue){
                minValue = inputArray[i];
            }
        }
        return minValue;
    }

}
