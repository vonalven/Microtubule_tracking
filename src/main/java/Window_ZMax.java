import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

public class Window_ZMax implements PlugIn {

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
        dlg.addNumericField("n stack proj", 1, 0);
        dlg.showDialog();
        if (dlg.wasCanceled()) return;
        double nStackProj = dlg.getNextNumber();

        nt = imp.getStackSize();
        int width = imp.getWidth();
        int height = imp.getHeight();

        int nOutStacks = (int)(nt-nStackProj);
        ImagePlus out = IJ.createHyperStack("Correction", width, height, 1, 1, nOutStacks, 16);
        out.getProcessor().resetMinAndMax();
        out.show();

        for (int t = 1; t <= nOutStacks; t++) {

            double maxIntensities[][] = new double[width][height]; // 0 by default

            for (int tt = 0; tt<nStackProj; tt++){
                imp.setPosition(t+tt);
                ImageProcessor ip = imp.getProcessor();

                for (int x = 0; x < width; x++) {
                    for (int y = 0; y < height; y++) {
                        maxIntensities[x][y] = Math.max(ip.getPixelValue(x, y),maxIntensities[x][y]);
                    }
                }
            }

            out.setPositionWithoutUpdate(1, 1, t);
            ImageProcessor op = out.getProcessor();
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    double v = (maxIntensities[x][y]);
                    op.putPixelValue(x, y, v);
                }
            }
        }

        if(true){
            return;
        }

    }
}
