/* *****************************************************************************
 *  Name: SeamCarver
 *  Description: SeamCarver finds the path of lowest energy pixels on a picture
 *  and removes the pixels along the path.
 **************************************************************************** */

import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

public class SeamCarver {
    private static final double BORDER_ENERGY = 1000.0;
    private Picture picture;
    private double[][] energyArr;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException();
        this.picture = new Picture(picture);
        energyArr = new double[picture.height()][picture.width()];
        for (double[] row : energyArr) {
            Arrays.fill(row, -1.0);
        }
    }

    // current picture
    public Picture picture() {
        return new Picture(picture);
    }

    // width of current picture
    public int width() {
        return energyArr[0].length;
    }

    // height of current picture
    public int height() {
        return energyArr.length;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        areValidIndices(x, y);
        if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1) return BORDER_ENERGY;
        if (energyArr[y][x] == -1.0)
            energyArr[y][x] = Math.sqrt(xGradient(x, y) + yGradient(x, y));
        return energyArr[y][x];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        energyArr = transpose(energyArr);
        transposePic();
        int[] seam = findVerticalSeam();
        energyArr = transpose(energyArr);
        transposePic();
        return seam;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        int[][] nodeTo = new int[height()][width()];
        double[][] distTo = new double[height()][width()];
        int[] seam = new int[height()];
        for (double[] distRow : distTo) {
            Arrays.fill(distRow, Double.POSITIVE_INFINITY);
        }
        Arrays.fill(distTo[0], 0.0);
        for (int y = 0; y < height() - 1; y++) {
            for (int x = 0; x < width(); x++) {
                for (int i = -1; i < 2; i++) {
                    double engy = energy(x, y);
                    if (!(x + i < 0 || x + i >= width())
                            && distTo[y + 1][x + i] > distTo[y][x] + engy) {
                        distTo[y + 1][x + i] = distTo[y][x] + engy;
                        nodeTo[y + 1][x + i] = x;
                    }

                }
            }
        }
        // find x value of pixel in bottom row with lowest distance
        int index = minIndex(distTo[height() - 1]);
        // loop through nodeTo and add path to seam[];
        seam[height() - 1] = index;
        for (int i = height() - 2; i >= 0; i--) {
            seam[i] = nodeTo[i + 1][index];
            index = seam[i];
        }

        return seam;
    }

    // finds the index of the double array element with the lowest value
    private int minIndex(double[] arr) {
        double min = Double.POSITIVE_INFINITY;
        int index = -1;
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] < min) {
                min = arr[i];
                index = i;
            }
        }
        return index;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        energyArr = transpose(energyArr);
        transposePic();
        removeVerticalSeam(seam);
        energyArr = transpose(energyArr);
        transposePic();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        isValidSeam(seam);
        Picture carvedPicture = new Picture(width() - 1, height());
        double[][] carvedEnergyArr = new double[height()][width() - 1];
        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < seam[y]; x++) {
                carvedPicture.set(x, y, picture.get(x, y));
                carvedEnergyArr[y][x] = energyArr[y][x];
            }
            if (seam[y] > 0)
                carvedEnergyArr[y][seam[y] - 1] = -1.0;
            for (int x = seam[y] + 1; x < width(); x++) {
                carvedPicture.set(x - 1, y, picture.get(x, y));
                carvedEnergyArr[y][x - 1] = energyArr[y][x];
            }
            if (seam[y] + 1 < width())
                carvedEnergyArr[y][seam[y]] = -1.0;
        }
        picture = carvedPicture;
        energyArr = carvedEnergyArr;
    }

    private void isValidSeam(int[] seam) {
        if (seam == null) throw new NullPointerException("Seam is null");
        if (seam.length != height())
            throw new IllegalArgumentException("seam length not equal to picture height");
    }

    private void areValidIndices(int x, int y) {
        if (x < 0 || x >= width() || y < 0 || y >= height())
            throw new IllegalArgumentException("Invalid indices");
    }

    private int xGradient(int x, int y) {
        int[] rgbLeft = extractRGBComp(picture.getRGB(x - 1, y));
        int[] rgbRight = extractRGBComp(picture.getRGB(x + 1, y));
        return dualGradient(rgbLeft, rgbRight);
    }

    private int yGradient(int x, int y) {
        int[] rgbUp = extractRGBComp(picture.getRGB(x, y - 1));
        int[] rgbDown = extractRGBComp(picture.getRGB(x, y + 1));
        return dualGradient(rgbUp, rgbDown);
    }

    private int dualGradient(int[] rgb1, int[] rgb2) {
        int gradient = 0;
        for (int i = 0; i < rgb1.length; i++) {
            gradient += Math.pow(rgb1[i] - rgb2[i], 2);
        }
        return gradient;
    }

    private int[] extractRGBComp(int rgb) {
        int[] rgbComps = new int[3];
        rgbComps[0] = (rgb >> 16) & 0xFF;
        rgbComps[1] = (rgb >> 8) & 0xFF;
        rgbComps[2] = rgb & 0xFF;
        return rgbComps;
    }


    private double[][] transpose(double[][] matrix) {
        double[][] T = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                T[j][i] = matrix[i][j];
            }
        }
        return T;
    }

    private void transposePic() {
        Picture transposedPic = new Picture(picture.height(), picture.width());
        for (int i = 0; i < picture.height(); i++) {
            for (int j = 0; j < picture.width(); j++) {
                transposedPic.setRGB(i, j, picture.getRGB(j, i));
            }
        }
        picture = transposedPic;
    }

    public static void main(String[] args) {
        Picture picture = new Picture(args[0]);
        SeamCarver carver = new SeamCarver(picture);
        int[] verticalSeam = carver.findVerticalSeam();
        for (int x : verticalSeam)
            StdOut.print(x + " ");
        StdOut.println("}");
    }
}
