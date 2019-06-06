import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class PuntaCometa {

    private int nt;

    public void run(int deltaFrame,int wNeighbors) {

        ImagePlus imp = IJ.getImage();
        ImagePlus imp2 = imp.duplicate();

        if (imp == null) {
            IJ.error("No open sequence of images.");
            return;
        }
        if (imp.getNFrames() < 1) {
            IJ.error("Requires a image with at least 2 frames.");
            return;
        }

        nt = imp.getStackSize();
        int width = imp.getWidth();
        int height = imp.getHeight();

        int nOutStacks = (nt-deltaFrame);
        ImagePlus out = IJ.createHyperStack("CometPoint2", width, height, 1, 1, nOutStacks, 16);
        out.getProcessor().resetMinAndMax();
        out.show();

        for (int t = 1; t <= nOutStacks; t++) {

            imp.setPosition(t);
            ImageProcessor ip1 = imp.getProcessor();
            double max1 = ip1.getMax();

            imp2.setPosition(t+deltaFrame);
            ImageProcessor ip2 = imp2.getProcessor();
            double max2 = ip2.getMax();

            double normMax = (max1+max2)/2;

            out.setPositionWithoutUpdate(1, 1, t);
            ImageProcessor op = out.getProcessor();

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    double neigMean = getNeighborsMean(ip1,x,y,wNeighbors);
                    int v = (int)(ip1.getPixelValue(x,y)*neigMean/normMax);
                    if(v<0){
                        v = 0;
                    }
                    op.putPixelValue(x,y,v);
                }
            }
        }

        ij.IJ.run("Enhance Contrast", "saturated=0.35");

    }



    private double[][] getNeighbors(ImageProcessor ip, int px, int py,int w){
        double[][] neighbors = new double[w*2+1][w*2+1];
        for (int x = px-w; x < px+w; x++) {
            for (int y = py-w; y < py+w; y++) {
                neighbors[x-px+w][y-py+w] = ip.getPixelValue(x,y);
            }
        }
        return neighbors;
    }

    private double getMean(double[][] inputArray){
        int wx = inputArray.length;
        int wy = inputArray[0].length;

        double mean = 0;

        for (int i = 0; i < wx; i++) {
            for (int j = 0; j < wy; j++) {
                mean+=inputArray[i][j];
            }
        }

        return mean/(wx*wy);
    }

    private double getNeighborsMean(ImageProcessor ip, int px, int py,int w){
        double[][] neighboars = getNeighbors(ip, px, py,w);
        return getMean(neighboars);
    }


}
