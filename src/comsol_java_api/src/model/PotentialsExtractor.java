package model;

import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import com.comsol.util.exceptions.FlException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;

/**
 * model.PotentialsExtractor
 *
 */
public class PotentialsExtractor {

    // INSTANCE VARIABLES
    private Model model; // model

    public IdentifierManager im = new IdentifierManager(); // top level identifier manager
    public static IdentifierManager ve_im = new IdentifierManager(); // top level identifier manager

    // directory structure
    private String root;

    // CONSTRUCTORS
    /**
     * Default constructor (minimum of 2 arguments)
     *
     * @param model       com.comsol.model.Model object is REQUIRED
     * @param projectRoot the root directory of the project (might remove if
     *                    unnecessary)
     */
    PotentialsExtractor(Model model, String projectRoot) {
        this.model = model;
        this.root = projectRoot;
    }

    // ACCESSOR/MUTATOR METHODS
    /**
     * @return the model
     */
    public Model getModel() {
        return model;
    }

    /**
     * @return the root of the project (String path)
     */
    public String getRoot() {
        return root;
    }

    /**
     * @param root set the project root (String path)
     */
    public void setRoot(String root) {
        this.root = root;
    }

    // OTHER METHODS
    /**
     * call method on im (IdentifierManager)... see class for details
     */
    public String next(String key) {
        return this.im.next(key);
    }

    /**
     * call method on im (IdentifierManager)... see class for details
     */
    public String next(String key, String pseudonym) {
        return this.im.next(key, pseudonym);
    }

    /**
     * @return success indicator
     */
    public static double[] extractPotentials(Model model, String coords_path) throws IOException {
        // Load coordinates (x,y,z) from file in form: top line is number of rows of
        // coords (int)
        // coordinates[0][i] = [x] in micron, (double)
        // coordinates[1][i] = [y] in micron, (double)
        // coordinates[2][i] = [z] in micron (double)

        // read in coords for axon segments as defined and saved to file in Python
        double[][] coordinatesLoaded;
        coordinatesLoaded = readCoords(coords_path);

        // transpose saved coordinates (we like to save (x,y,z) as column vectors, but
        // COMSOL wants as rows)
        double[][] coordinates;
        coordinates = transposeMatrix(coordinatesLoaded);

        // get Ve from COMSOL
        String id = ve_im.next("interp");
        // Check each dataset for V
        String[] tags = model.result().dataset().tags();
        for (int tagno = 0; tagno < tags.length; tagno++) {
            try {
                model.result().numerical().create(id, "Interp");
                model.result().numerical(id).set("expr", "V");
                model.result().numerical(id).set("data", tags[tagno]);
                model.result().numerical(id).setInterpolationCoordinates(coordinates);

                double[][][] ve_pre = model.result().numerical(id).getData();
                int len = ve_pre[0][0].length; // number of coordinates

                double[] ve = new double[len];
                for (int i = 0; i < len; i++) {
                    ve[i] = ve_pre[0][0][i];
                }
                //System.out.println("Successfully extracted V from dataset "+tags[tagno]+".");
                return ve;
            } catch (Exception e) {
                // System.out.println("Could not extract V for dataset "+tags[tagno]+".");
                model.result().numerical().remove(id);
                if (tagno == tags.length - 1) {
                    throw e;
                }
            }
        }
        return new double[0];
    }

    // https://stackoverflow.com/questions/15449711/transpose-double-matrix-with-a-java-function
    public static double[][] transposeMatrix(double[][] m) {
        // pre-allocated array of doubles for transposed matrix
        double[][] temp = new double[m[0].length][m.length];

        for (int i = 0; i < m.length; i++) for (int j = 0; j < m[0].length; j++) temp[j][i] =
            m[i][j];
        return temp;
    }

    private static boolean writeVe(double[] ve, String ve_path) throws IOException {
        PrintWriter printWriter = new PrintWriter(ve_path);
        int len = ve.length; // number of coordinates

        // write to file: number of coordintates top line,
        // then one Ve value for each coordinate (x,y,z) for subsequent lines
        printWriter.println(len);
        for (int i = 0; i < len; i++) {
            printWriter.println(ve[i]);
        }
        printWriter.close(); // close printWriter
        return true;
    }

    public static double[][] readCoords(String coords_path) throws FileNotFoundException {
        File f = new File(coords_path);
        Scanner scan = new Scanner(f);

        String thisLine = null;
        try {
            // save rows (number of coords) at top line... so number of lines in file is
            // (number of coords +1)
            String rows = scan.nextLine();
            int n_rows = Integer.parseInt(rows.trim());

            // pre-allocated array of doubles for coords in file (3 columns by default for
            // (x,y,z)
            double[][] coords = new double[n_rows][3];
            int row_ind = 0;

            // while there are more lines to scan
            while (scan.hasNextLine()) {
                thisLine = scan.nextLine();
                String[] parts = thisLine.split("\\s+");
                for (int i = 0; i < parts.length; i++) {
                    coords[row_ind][i] = Double.parseDouble(parts[i]);
                }
                row_ind++;
            }

            scan.close();

            if (n_rows != row_ind) {
                throw new Exception(
                    "Number of coordinates (rows) in coords file does not match header in file: " +
                    coords_path
                );
            }

            return coords;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

    /**
     * Master procedure to run!
     *
     * @param args
     */
    public static void main(String[] args) throws InterruptedException {
        // Start COMSOL Instance
        try {
            ModelUtil.connect("localhost", 2036);
        } catch (FlException e) {
            ModelUtil.connect("localhost", 2037);
        }

        TimeUnit.SECONDS.sleep(5);
        System.out.println("Connected to COMSOL Server");
        ModelUtil.initStandalone(false);
        // ModelUtil.showProgress(null);

        // Take projectPath input to PotentialsExtractor and assign to string.
        String mphFile = args[0];

        // Load RUN configuration data
        String fibersPath = args[1]; // Take runPath input to PotentialsExtractor and assign to string

        String outPath = args[2];

        System.out.println("loading mph file: " + mphFile);

        // Define PotentialsExtractor class instance for model and projectPath
        Model model = null;
        try {
            model = ModelUtil.load("Model", mphFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        // PotentialsExtractor mw = new PotentialsExtractor(model, projectPath);
        File dir = new File(fibersPath);
        File[] directoryListing = dir.listFiles();
        if (directoryListing != null) {
            for (File child : directoryListing) {
                String coords_file = child.getAbsolutePath();
                if (coords_file.endsWith(".dat")) {
                    //System.out.println("loading coords file: " + coords_file);
                    double[] ve = null;
                    try {
                        ve = extractPotentials(model, coords_file);
                    } catch (Exception e) {
                        System.out.println("Failed to extract potentials.");
                        System.out.println(e);
                        continue;
                    }

                    String ve_file_path = String.join(
                        "/",
                        new String[] { outPath + "/" + child.getName() }
                    );
                    try {
                        writeVe(ve, ve_file_path);
                        //System.out.println("Saved ve to: " + ve_file_path);
                    } catch (IOException e) {
                        System.out.println("Failed to extract potentials.");
                        e.printStackTrace();
                    }
                }
            }
        }
        ModelUtil.remove(model.tag());
        ModelUtil.disconnect();
        System.out.println("Disconnected from COMSOL Server");
        System.exit(0);
    }
}
