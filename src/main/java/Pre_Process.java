import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class Pre_Process implements PlugIn {

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

        ij.IJ.run("Duplicate...", "title=pre_process duplicate");
        ij.IJ.run("Subtract Background...", "rolling=20 sliding stack");
        ij.IJ.run("Enhance Contrast", "saturated=0.35");
        ij.IJ.run("Window ZMax", "n=1");
        ij.IJ.run("Enhance Contrast", "saturated=0.35");
        ij.IJ.selectWindow("pre_process");
        ij.IJ.run("Close");
        ij.IJ.selectWindow("ZMaxWindow"); //TODO: rename this window
        ij.IJ.run("Sigma Filter Plus", "radius=1 use=2 minimum=0.2 outlier stack");
        ij.IJ.run("Subtract Background...", "rolling=50 sliding stack");
        ij.IJ.run("Enhance Contrast", "saturated=0.35");
        ij.IJ.run("Punta Cometa2", "delta frame=7 neighbors=3");
        ij.IJ.selectWindow("ZMaxWindow");
        ij.IJ.run("Close");
        ij.IJ.selectWindow("CometPoint2");
    }
}