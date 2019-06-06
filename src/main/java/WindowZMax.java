import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

public class WindowZMax {

    private int nt;

    public void run(double nStackProj) {

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
        int width = imp.getWidth();
        int height = imp.getHeight();

        int nOutStacks = (int)(nt-nStackProj);
        ImagePlus out = IJ.createHyperStack("ZMaxWindow", width, height, 1, 1, nOutStacks, 16);
        out.getProcessor().resetMinAndMax();
        out.show();

        for (int t = 1; t <= nOutStacks; t++) {

            double[][] maxIntensities = new double[width][height]; // 0 by default

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

    }
}
