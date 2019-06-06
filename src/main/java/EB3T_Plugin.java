import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class EB3T_Plugin implements PlugIn {

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

        GenericDialog dlg = new GenericDialog("Particle Tracking");
        dlg.addMessage("Welcome to the EB3Tracker");
        dlg.addCheckbox("Pre Process image", false);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        boolean preProcess = dlg.getNextBoolean();

        if(preProcess) {
            ij.IJ.run("EB3T Pre Process");
        }
        ij.IJ.run("EB3T Particle Tracking");

    }

    public static void main(String[] args) {

    }
}