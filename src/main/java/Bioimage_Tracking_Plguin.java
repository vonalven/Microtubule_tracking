import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class Bioimage_Tracking_Plguin implements PlugIn {

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
        ij.IJ.run("Particle Tracking", "threshold=0.4 radius=2.5 preview max=10 min=7");

    }

    public static void main(String[] args) {

    }
}