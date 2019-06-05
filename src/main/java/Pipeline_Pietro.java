import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class Pipeline_Pietro implements PlugIn {

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

        ij.IJ.run("Pre Process");
        ij.IJ.run("Auto SandT");

    }
}