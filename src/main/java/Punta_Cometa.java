import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class Punta_Cometa implements PlugIn {

    private int nt;

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

        GenericDialog dlg = new GenericDialog("TEST");
        dlg.addNumericField("threshold", 20, 3);
        dlg.addNumericField("delta frame", 2, 0);
        dlg.addNumericField("w size", 2, 0);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        double threshold = dlg.getNextNumber();
        int deltaFrame = (int)dlg.getNextNumber();
        int wsize = (int)dlg.getNextNumber();

        nt = imp.getStackSize();
        int width = imp.getWidth();
        int height = imp.getHeight();

        int nOutStacks = (int)(nt-deltaFrame);
        ImagePlus out = IJ.createHyperStack("CometPoint", width, height, 1, 1, nOutStacks, 16);
        out.getProcessor().resetMinAndMax();
        out.show();

        double meanBefore = 0;
        double meanAfter  = 0;

        for (int t = 1; t <= nOutStacks; t++) {

            imp.setPosition(t);
            ImageProcessor ip1 = imp.getProcessor();

            ImagePlus imp2 = imp.duplicate();

            imp2.setPosition(t+deltaFrame);
            ImageProcessor ip2 = imp2.getProcessor();

            out.setPositionWithoutUpdate(1, 1, t);
            ImageProcessor op = out.getProcessor();


            for (int x = wsize; x < width-wsize; x++) {
                for (int y = wsize; y < height-wsize; y++) {

                    meanBefore = getMean(getNeighbors(ip1,x,y,wsize));
                    meanAfter = getMean(getNeighbors(ip2,x,y,wsize));

                    if(Math.abs(meanAfter-meanBefore)<threshold){
                        op.putPixelValue(x, y, Math.round(meanAfter-meanBefore)/2);
                    }else{
                        op.putPixelValue(x, y, ip1.getPixelValue(x,y));
                    }

                }
            }
        }

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

}
