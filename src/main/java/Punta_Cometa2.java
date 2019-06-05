import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class Punta_Cometa2 implements PlugIn {

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
        dlg.addNumericField("delta frame", 5, 0);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        int deltaFrame = (int)dlg.getNextNumber();

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

            ImagePlus imp2 = imp.duplicate();

            imp2.setPosition(t+deltaFrame);
            ImageProcessor ip2 = imp2.getProcessor();

            out.setPositionWithoutUpdate(1, 1, t);
            ImageProcessor op = out.getProcessor();

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    int v = (int)(ip1.getPixelValue(x,y)-ip2.getPixelValue(x,y));
                    if(v<0){
                        v = 0;
                    }
                    op.putPixelValue(x,y,v);
                }
            }
        }

    }

}
