/*
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
*/

package threedmodel;

import com.comsol.model.*;
import com.comsol.model.util.ModelUtil;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * model.ModelWrapper
 * Master high-level class for managing a model, its metadata, and various critical operations such as creating parts,
 * assigning physics, and extracting potentials. This class houses the "meaty" operations of actually interacting with
 * the model object when creating parts in the static class model.Parts.
 */
public class ModelWrapper {

    // UNION PSEUDONYM CONSTANTS
    public static final String ALL_NERVE_PARTS_UNION = "allNervePartsUnion";
    public static final String ENDO_UNION = "endoUnion";
    public static final String PERI_UNION = "periUnion";

    // if these change, also need to change in createEnvironmentPartInstance
    public static final String DISTAL_MEDIUM = "DistalMedium";
    public static final String PROXIMAL_MEDIUM = "ProximalMedium";

    public static final String[] ALL_UNIONS = new String[] {
        ModelWrapper.ENDO_UNION,
        ModelWrapper.ALL_NERVE_PARTS_UNION,
        ModelWrapper.PERI_UNION,
    };
    // associated union contributors for above constants
    private final HashMap<String, ArrayList<String>> unionContributors = new HashMap<>();

    // INSTANCE VARIABLES
    private final Model model; // model

    public IdentifierManager im = new IdentifierManager(); // top level identifier manager

    private final HashMap<String, IdentifierManager> partPrimitiveIMs = new HashMap<>(); // for managing parts within COMSOL

    // directory structure
    private String root;
    private String dest;

    // CONSTRUCTORS

    /**
     * Default constructor (minimum of 2 arguments)
     *
     * @param model       com.comsol.model.Model object is REQUIRED
     * @param projectRoot the root directory of the project (might remove if unnecessary)
     */
    ModelWrapper(Model model, String projectRoot) {
        this.model = model;
        this.root = projectRoot;
        this.initUnionContributors();
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
     * call method on im (IdentifierManager)... see class for details
     */
    public String get(String psuedonym) {
        return this.im.get(psuedonym);
    }

    /**
     * @param partPrimitiveLabel the name of the part primitive (i.e. "TubeCuff_Primitive")
     * @return the associated IdentifierManager, for correct intra-part indexing
     */
    public IdentifierManager getPartPrimitiveIM(String partPrimitiveLabel) {
        return this.partPrimitiveIMs.get(partPrimitiveLabel);
    }

    /**
     * @param destination full path to save to
     * @return success indicator
     */
    public boolean save(String destination) {
        try {
            this.model.save(destination);
            return true;
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Convenience method for saving to relative directory (this.dest) wrt the project directory (root)
     *
     * @return success indicator
     */
    public boolean save() {
        if (this.dest != null) return save(
            String.join("/", new String[] { this.root, this.dest })
        ); else {
            System.out.println("Save directory not initialized");
            return false;
        }
    }

    /**
     * Create the required primitives for a given cuff json
     *
     * @param name json filename WITH extension (i.e. "LivaNova2000.json")
     */
    public void addCuffPartPrimitives(String name) {
        // extract data from json
        try {
            JSONObject cuffData = JSONio.read(
                String.join("/", new String[] { this.root, "config", "system", "cuffs", name })
            );

            // get the id for the next "par" (i.e., parameters section), and give it a name from the JSON file name
            String id = this.next("par", name);
            model.param().group().create(id);
            model.param(id).label(name.split("\\.")[0] + " Parameters");

            // loop through all parameters in file, and set in parameters
            for (Object item : (JSONArray) cuffData.get("params")) {
                JSONObject itemObject = (JSONObject) item;
                model
                    .param(id)
                    .set(
                        (String) itemObject.get("name"),
                        (String) itemObject.get("expression"),
                        (String) itemObject.get("description")
                    );
            }

            // for each required part primitive, create it (if not already existing)
            for (Object item : (JSONArray) cuffData.get("instances")) {
                JSONObject itemObject = (JSONObject) item;
                String partPrimitiveName = (String) itemObject.get("type"); // quick cast to String

                // create the part primitive if it has not already been created
                if (!this.im.hasPseudonym(partPrimitiveName)) {
                    // get next available (TOP LEVEL) "part" id
                    String partID = this.im.next("part", partPrimitiveName);
                    try {
                        // TRY to create the part primitive (catch error if no existing implementation)
                        IdentifierManager partPrimitiveIM = Part.createCuffPartPrimitive(
                            partID,
                            partPrimitiveName,
                            this,
                            (JSONArray) cuffData.get("local_params")
                        );

                        // add the returned id manager to the HashMap of IMs with the partName as its key
                        this.partPrimitiveIMs.put(partPrimitiveName, partPrimitiveIM);
                    } catch (IllegalArgumentException e) {
                        e.printStackTrace();
                        return;
                    }
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
     * Instantiate required primitives for given cuff
     * NOTE: addCuffPartPrimitives() MUST be called first or there will be no primitives to instantiate
     *
     * @param name same formatting as in addCuffPartPrimitives()
     */
    public void addCuffPartInstances(String name, String pcsfolder) {
        // extract data from json (name is something like Enteromedics.json)
        try {
            JSONObject cuffData = JSONio.read(
                String.join("/", new String[] { this.root, "config", "system", "cuffs", name })
            );

            // loop through all part instances
            for (Object item : (JSONArray) cuffData.get("instances")) {
                JSONObject itemObject = (JSONObject) item;

                String instanceLabel = (String) itemObject.get("label");
                String instanceID = this.im.next("pi", instanceLabel);
                String type = (String) itemObject.get("type");
                Part.createCuffPartInstance(
                    instanceID,
                    instanceLabel,
                    type,
                    this,
                    itemObject,
                    pcsfolder
                );
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
     * Create materials necessary for fascicles, nerve, surrounding media, etc.
     */
    public void addMaterialDefinitions(
        ArrayList<String> materials,
        JSONObject modelData,
        ModelParamGroup materialParams
    ) {
        try {
            // load system defined materials JSON into memory
            JSONObject materialsData = JSONio.read(
                String.join("/", new String[] { this.root, "config", "system", "materials.json" })
            );

            // add material definition for those materials that are needed in the instantiated parts
            for (String function : materials) {
                if (!this.im.hasPseudonym(function)) {
                    String materialID = this.im.next("mat", function);
                    Part.defineMaterial(
                        materialID,
                        function,
                        modelData,
                        materialsData,
                        this,
                        materialParams
                    );
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    /**
     * Call only from initializer!
     * Initialize the ArrayLists in the unionContributors HashMap
     */
    public void initUnionContributors() {
        for (String unionLabel : ModelWrapper.ALL_UNIONS) {
            this.unionContributors.put(unionLabel, new ArrayList<>());
        }
    }

    /**
     * Add string id for COMSOL element to the listed unions (which have not be "created" in COMSOL yet)
     *
     * @param contributor the string id to add (use actual id, not pseudonym)
     * @param unions      which unions to add it to  (use static pseudonym constants at top of class)
     */
    public void contributeToUnions(String contributor, String[] unions) {
        for (String union : unions) {
            this.unionContributors.get(union).add(contributor);
        }
    }

    /**
     * @param union which union to from of which to get the contributors
     * @return String array of the COMSOL id's of contributors (likely ext# or csel#)
     */
    public String[] getUnionContributors(String union) {
        if (!this.unionContributors.containsKey(union)) throw new IllegalArgumentException(
            "No such union: " + union
        );
        return this.unionContributors.get(union).toArray(new String[0]);
    }

    /**
     * Actually create the unions by looping through all defined ArrayLists and adding contents to a new union.
     * Will not create a union of no elements in associated ArrayList (i.e. no Peri union if only contact impedance)
     */
    public void createUnions() {
        for (String union : ModelWrapper.ALL_UNIONS) {
            String[] contributors = this.getUnionContributors(union);

            if (contributors.length > 0) {
                GeomFeature uni = model
                    .component("comp1")
                    .geom("geom1")
                    .create(im.next("uni", union), "Union");
                uni.set("keep", true);
                uni.selection("input").set(contributors);
                uni.label(union);

                String unionCselLabel = union + "Csel";
                GeomObjectSelectionFeature csel = model
                    .component("comp1")
                    .geom("geom1")
                    .selection()
                    .create(im.next("csel", unionCselLabel), "CumulativeSelection");
                csel.label(unionCselLabel);

                uni.set("contributeto", im.get(unionCselLabel));
            }
        }
    }

    /**
     * Master procedure to run!
     */

    public static void main(String[] args) throws InterruptedException, IOException {
        // Start COMSOL Instance
        try {
            ModelUtil.connect();
        } catch (Exception e) {
            System.out.println("COMSOL connection failed, retrying with port 2036...");
            try {
                ModelUtil.connect("localhost", 2036);
            } catch (Exception e2) {
                System.out.println("COMSOL connection failed, retrying with specific port...");
                try {
                    ModelUtil.connect("localhost", 2037);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
            }
        }

        TimeUnit.SECONDS.sleep(5);
        ModelUtil.initStandalone(false);
        //        ModelUtil.showProgress(null); // if you want to see COMSOL progress (as it makes all geometry, runs, etc.)

        // Take projectPath input to ModelWrapper and assign to string.
        String projectPath = args[0];

        JSONObject config3D = null;
        try {
            config3D = JSONio.read(args[1]);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        String runtype = args[2];

        long wait_hours = 24;

        if (runtype.equals("doubledeform")) {
            LicenseCheckout(wait_hours, "STRUCTURALMECHANICS");
        }

        String licenseType = "COMSOL";
        LicenseCheckout(wait_hours, licenseType);

        // Take runPath input to ModelWrapper and assign to string
        Model model = null;
        ModelWrapper mw;

        // Load MODEL configuration data

        JSONObject modelData = null;
        try {
            assert config3D != null;
            modelData = JSONio.read(config3D.getString("project_path") + "/model.json");
        } catch (FileNotFoundException e) {
            System.out.println("Failed to read MODEL config data.");
            e.printStackTrace();
        }
        switch (runtype) {
            case "geometry":
                model = buildGeometry(config3D, modelData, projectPath);
                break;
            //Run deformation algorithm
            case "doubledeform":
                model = runDeformation(config3D);
                break;
            //Generate electric currents simulation in COMSOL
            case "model":
                model = solveModels(config3D, modelData, projectPath);
                break;
            //Generate electric currents simulation in COMSOL
            case "fibers":
                model = getFibers(config3D);
                break;
            default:
                System.out.println("INVALID CHOICE FOR RUN TYPE");
                break;
        }
        assert model != null;
        ModelUtil.remove(model.tag());

        System.exit(0);
    }

    private static Model getFibers(JSONObject config3D) throws IOException {
        Model model;
        System.out.println("Generating Fibers");
        ModelUtil.showProgress(null);
        String model_tag = ModelUtil.uniquetag("Model");

        String outpath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("comsol");
        String geomFile = outpath + "/0.mph";
        String fiberFile = outpath + "/fibermodel.mph";

        model = ModelUtil.load(model_tag, geomFile);

        String fiberPath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("fibers");
        // START PRE MESH
        System.out.println("Setting up fiber model...");

        //Generate fiberset
        model.result().create("pg1", "PlotGroup3D");
        model.result("pg1").label("FiberSet");
        model.result("pg1").set("data", "dset1");
        model.result("pg1").selection().geom("geom1", 3);
        model
            .result("pg1")
            .selection()
            .set(model.component("comp1").physics("cc").selection().entities(3));
        model.result("pg1").create("str1", "Streamline");
        model.result("pg1").feature("str1").set("selnumber", 250);
        model.result("pg1").feature("str1").set("expr", new String[] { "cc.vX", "cc.vY", "cc.vZ" });
        model.result("pg1").feature("str1").set("maxlen", Double.POSITIVE_INFINITY);
        model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
        model.result("pg1").feature("str1").set("posmethod", "selection");
        model.result("pg1").feature("str1").selection().named("box2");
        model.result("pg1").feature("str1").set("linetype", "line");
        model.result("pg1").feature("str1").set("data", "dset1");
        model.result("pg1").feature("str1").set("resolution", "normal");
        model.result("pg1").feature("str1").set("smooth", "material");
        model.result("pg1").feature("str1").set("recover", false);
        model.result("pg1").feature("str1").set("inttol", "0.00001");
        model.result("pg1").feature("str1").set("maxsteps", 50000);
        model.result("pg1").feature("str1").set("stattol", "0.0001");
        model.result().export().create("plot1", "pg1", "str1", "Plot");
        model.result().export("plot1").set("filename", fiberPath + "/fibers.txt");
        model.result().export("plot1").set("struct", "spreadsheet");
        model.result().export("plot1").set("separator", ",");
        model.result("pg1").feature("str1").create("filt1", "Filter");
        model.result().export().create("plot2", "pg1", "str1", "Plot");
        model.result().export("plot2").set("compact", true);
        model.result().export("plot2").set("separator", ",");
        model.save(outpath + "/debug_fibermodel.mph");
        model.result("pg1").run();
        model.result().export("plot1").run();

        JSONObject zlocs = JSONio.read(fiberPath + "/zlocs.json");
        Iterator<String> zkeys = zlocs.keys();

        while (zkeys.hasNext()) {
            String key = zkeys.next();
            double microns = zlocs.getDouble(key);
            model
                .result("pg1")
                .feature("str1")
                .feature("filt1")
                .set("expr", "z>" + microns / 1000000);
            model.result().export("plot2").set("filename", fiberPath + "/" + key + ".txt");
            model.result("pg1").run();
            model.result().export("plot2").run();
        }
        model.save(fiberFile);
        return model;
    }

    private static Model solveModels(JSONObject config3D, JSONObject modelData, String projectPath)
        throws IOException {
        ModelWrapper mw;
        Model model;
        String outpath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("comsol");
        String geomFile = outpath + "/model.mph";

        // START PRE MESH
        System.out.println("Setting up model...");

        // Define model object
        model = ModelUtil.createUnique("Model");
        // Add component node 1
        model.component().create("comp1", true);
        // Add 3D geom to component node 1
        model.component("comp1").geom().create("geom1", 3);
        // geometry shape order
        assert modelData != null;
        String sorder = modelData.getJSONObject("solver").getString("sorder");
        model.component("comp1").sorder(sorder);
        // Set default length units to micron
        model.component("comp1").geom("geom1").lengthUnit("\u00b5m");
        // Add materials node to component node 1
        model.component("comp1").physics().create("ec", "ConductiveMedia", "geom1");
        // and mesh node to component node 1
        model.component("comp1").mesh().create("mesh1");

        // Define ModelWrapper class instance for model and projectPath
        mw = new ModelWrapper(model, projectPath);

        String meshpath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("mesh") +
            "/mesh.nas";
        String pcspath = outpath + "/pcs";
        // Add materials
        System.out.println("generating material params");
        // add MATERIAL DEFINITIONS
        String materialParamsLabel = "Material Parameters";
        ModelParamGroup materialParams = model.param().group().create(materialParamsLabel);
        materialParams.label(materialParamsLabel);

        ArrayList<String> bio_materials = new ArrayList<>(
            Arrays.asList("medium", "perineurium", "endoneurium", "epineurium")
        );

        JSONObject cuffObject = (JSONObject) modelData.get("cuff");
        String cuff = cuffObject.getString("preset");
        JSONObject cuffData = JSONio.read(
            String.join("/", new String[] { mw.root, "config", "system", "cuffs", cuff })
        );

        ArrayList<String> cuff_materials = new ArrayList<>();
        // loop through all part instances
        for (Object item : (JSONArray) cuffData.get("instances")) {
            JSONObject itemObject = (JSONObject) item;
            for (Object function : itemObject.getJSONArray("materials")) {
                JSONObject functionObject = (JSONObject) function;
                cuff_materials.add(functionObject.getString("info"));
            }
        }

        mw.addMaterialDefinitions(cuff_materials, modelData, materialParams);

        mw.addMaterialDefinitions(bio_materials, modelData, materialParams);

        System.out.println("importing meshes");
        System.out.println(meshpath);

        //define numpcs as the number of pcs files in the pcs folder
        int numpcs = 0;
        File pcsdir = new File(pcspath);
        File[] pcsfiles = pcsdir.listFiles();
        assert pcsfiles != null;
        for (File file : pcsfiles) {
            if (file.getName().endsWith(".dat")) {
                numpcs++;
            }
        }
        double[][] pcsPoints = new double[numpcs][3];
        for (int i = 0; i < numpcs; i++) {
            double[] pt = new double[3];
            try (
                BufferedReader br = new BufferedReader(
                    new FileReader(pcspath + "/pcs" + (i + 1) + ".dat")
                )
            ) {
                for (int j = 0; j < 3; j++) {
                    String line = br.readLine();
                    pt[j] = Double.parseDouble(line);
                }
            }
            pcsPoints[i] = pt;
        }
        try {
            //import mesh and add necessary steps
            model.component("comp1").mesh("mesh1").create("imp1", "Import");
            model.component("comp1").mesh("mesh1").feature("imp1").set("filename", meshpath);
            model.component("comp1").mesh("mesh1").feature("imp1").set("facepartition", "minimal");

            for (int i = 0; i < numpcs; i++) {
                model.component("comp1").mesh("mesh1").create("vtx" + (i + 1), "CreateVertices");
                model
                    .component("comp1")
                    .mesh("mesh1")
                    .feature("vtx" + (i + 1))
                    .set("vertexspec", "coord");
                model
                    .component("comp1")
                    .mesh("mesh1")
                    .feature("vtx" + (i + 1))
                    .set("x", pcsPoints[i][0]);
                model
                    .component("comp1")
                    .mesh("mesh1")
                    .feature("vtx" + (i + 1))
                    .set("y", pcsPoints[i][1]);
                model
                    .component("comp1")
                    .mesh("mesh1")
                    .feature("vtx" + (i + 1))
                    .set("z", pcsPoints[i][2]);
                model
                    .component("comp1")
                    .mesh("mesh1")
                    .feature("vtx" + (i + 1))
                    .set("relsnaptol", ".01");
            }

            model.component("comp1").mesh("mesh1").run("fin");
        } catch (Exception e) {
            System.out.println("An error occurred, saving model file for debugging");
            model.save(geomFile);
            e.printStackTrace();
        }

        System.out.println(
            "Mesh import complete, setting up physics and assigning boundary conditions"
        );

        //add curvilinear coordinates
        model.component("comp1").physics().create("cc", "CurvilinearCoordinates", "geom1");
        model.component("comp1").physics("cc").prop("Settings").set("CreateBasis", true);
        model.component("comp1").physics("cc").create("diff1", "DiffusionMethod", 3);
        model.component("comp1").physics("cc").feature("diff1").create("inl1", "Inlet", 2);
        model.component("comp1").physics("cc").feature("diff1").create("out1", "Outlet", 2);

        //Add current physics
        model.component("comp1").physics("ec").create("cucn2", "CurrentConservation", 3);
        model.component("comp1").physics("ec").feature("cucn2").set("coordinateSystem", "cc_cs");
        model
            .component("comp1")
            .physics("ec")
            .feature("cucn2")
            .set("ConstitutiveRelationJcE", "ElectricalConductivity");
        model.component("comp1").physics("ec").feature("cucn2").set("sigma_mat", "userdef");
        model
            .component("comp1")
            .physics("ec")
            .feature("cucn2")
            .set(
                "sigma",
                new String[] {
                    "sigma_endoneurium_z",
                    "0",
                    "0",
                    "0",
                    "sigma_endoneurium_x",
                    "0",
                    "0",
                    "0",
                    "sigma_endoneurium_y",
                }
            );
        model
            .component("comp1")
            .physics("ec")
            .feature("cucn2")
            .label("Curvilinear Current Conservation");
        model.component("comp1").physics("ec").create("gnd1", "Ground", 2);
        for (int i = 1; i < numpcs + 1; i++) {
            model.component("comp1").physics("ec").create("pcs" + i, "PointCurrentSource", 0);
        }
        model.component("comp1").selection().create("sel1", "Explicit");
        model.component("comp1").selection("sel1").geom("geom1", 3, 2, new String[] { "exterior" });
        model.component("comp1").selection("sel1").all();
        model.component("comp1").physics("ec").feature("gnd1").selection().named("sel1");
        model.save(geomFile);

        // Add material Links
        String[] sels = model.component("comp1").mesh("mesh1").feature("imp1").outputSelection();
        JSONObject meshmap = JSONio.read(outpath + "/stl/meshmap.json");
        JSONObject mats = meshmap.getJSONObject("materials");
        Iterator<String> keys = mats.keys();
        model.component("comp1").selection().create("com1", "Complement");
        List<String> sipsels = new ArrayList<>();
        List<String> pcsSels = new ArrayList<>();

        while (keys.hasNext()) {
            String matname = keys.next();
            JSONArray mat_domains = mats.getJSONArray(matname);
            for (int i = 0; i < mat_domains.length(); i++) {
                String id = mat_domains.getString(i);
                String linkLabel = matname + "_link_dom" + id;
                Material mat = model
                    .component("comp1")
                    .material()
                    .create(mw.im.next("matlnk", linkLabel), "Link");
                mat.label(linkLabel);
                mat.set("link", mw.im.get(matname));
                String this_nas = "imp1_nastransel_3_" + id;
                for (int j = 0; j < sels.length; j++) {
                    if (sels[j].equals(this_nas)) {
                        sipsels.add("nastransel" + (j + 1));
                        if (!(matname.equals("endoneurium"))) {
                            mat.selection().named("nastransel" + (j + 1));
                            if (matname.equals("conductor")) {
                                System.out.println("conductor added");
                                pcsSels.add("nastransel" + (j + 1));
                            }
                        } else {
                            model
                                .component("comp1")
                                .physics("cc")
                                .selection()
                                .named("nastransel" + (j + 1));
                            model
                                .component("comp1")
                                .physics("cc")
                                .feature("diff1")
                                .selection()
                                .named("nastransel" + (j + 1));
                            model
                                .component("comp1")
                                .physics("ec")
                                .feature("cucn2")
                                .selection()
                                .named("nastransel" + (j + 1));
                        }
                    }
                }
            }
        }

        double[] zpos = new double[numpcs];
        System.out.println(pcsSels);
        for (int i = 0; i < pcsSels.size(); i++) {
            model.component("comp1").selection().create("adj" + i, "Adjacent");
            model
                .component("comp1")
                .selection("adj" + i)
                .set("input", pcsSels.toArray(new String[pcsSels.size()])[i]);
            model.component("comp1").selection("adj" + i).set("outputdim", 0);
            model.component("comp1").selection("adj" + i).set("interior", true);
            model.component("comp1").selection("adj" + i).set("exterior", false);
            int[] ents = model.component("comp1").selection("adj" + i).entities(0);
            model.component("comp1").geom("geom1").measureFinal().selection().geom("geom1", 0);
            model.component("comp1").geom("geom1").measureFinal().selection().set(ents);
            double[] coordsel = model.component("comp1").geom("geom1").measureFinal().getVtxCoord();
            zpos[i] = coordsel[2];
        }
        System.out.println("success to here");
        //model.save("D:/debug.mph");

        //Assign all pcs
        //loop through each pcs point, and assign that pcs the Adj which is closest to it
        for (int i = 0; i < numpcs; i++) {
            Integer minAdj = null;
            //loop through all adj selections, equal count to pcs
            for (int adjnum = 0; adjnum < numpcs; adjnum++) {
                int[] ents = model.component("comp1").selection("adj" + adjnum).entities(0);
                model.component("comp1").measure().selection().geom("geom1", 0);
                //assert ents is length one
                assert ents.length == 1;
                model.component("comp1").measure().selection().set(ents);
                double[] adjCoord = model.component("comp1").measure().getVtxCoord();
                //calculate distance between adj and pcs
                double dist = Math.sqrt(
                    Math.pow(adjCoord[0] - pcsPoints[i][0], 2) +
                    Math.pow(adjCoord[1] - pcsPoints[i][1], 2) +
                    Math.pow(adjCoord[2] - pcsPoints[i][2], 2)
                );
                //if dist<100, assign pcs to adj
                if (dist < 100) {
                    //if minAdj is already set, error
                    if (minAdj != null) {
                        System.out.println("Error: multiple adjacent selections found for pcs");
                        System.exit(1);
                    }
                    minAdj = adjnum;
                }
            }
            model
                .component("comp1")
                .physics("ec")
                .feature("pcs" + (i + 1))
                .selection()
                .named("adj" + minAdj);
        }

        String matname = "medium";
        model
            .component("comp1")
            .selection("com1")
            .set("input", sipsels.toArray(new String[sipsels.size()]));
        String linkLabel = matname + "_link";
        Material mat = model
            .component("comp1")
            .material()
            .create(mw.im.next("matlnk", linkLabel), "Link");
        mat.label(linkLabel);
        mat.set("link", mw.im.get(matname));
        mat.selection().named("com1");

        JSONObject cons = meshmap.getJSONObject("contacts");
        String id_inout = cons.getString("i with medium");
        String sel_inout = "";
        for (int j = 0; j < sels.length; j++) {
            if (sels[j].equals("imp1_nastransel_2_" + id_inout)) {
                sel_inout = "nastransel" + (j + 1);
                System.out.println(sel_inout);
            }
        }
        model.component("comp1").selection().create("box1", "Box");
        model.component("comp1").selection("box1").set("entitydim", 2);
        model.component("comp1").selection("box1").set("inputent", "selections");
        model.component("comp1").selection("box1").set("input", new String[] { sel_inout });
        model.component("comp1").selection("box1").set("zmax", ".00125");
        model.component("comp1").selection("box1").set("condition", "inside");
        model.component("comp1").selection("box1").label("rostral");
        model.component("comp1").selection().create("box2", "Box");
        model.component("comp1").selection("box2").set("entitydim", 2);
        model.component("comp1").selection("box2").set("inputent", "selections");
        model.component("comp1").selection("box2").set("input", new String[] { sel_inout });
        model.component("comp1").selection("box2").set("zmin", ".00125");
        model.component("comp1").selection("box2").set("condition", "inside");
        model.component("comp1").selection("box2").label("caudal");

        model
            .component("comp1")
            .physics("cc")
            .feature("diff1")
            .feature("inl1")
            .selection()
            .named("box1");
        model
            .component("comp1")
            .physics("cc")
            .feature("diff1")
            .feature("out1")
            .selection()
            .named("box2");

        System.out.println("Saving pre-solution model");
        model.save(geomFile);

        //check that every pcs has one entity
        for (int pnum = 1; pnum < numpcs + 1; pnum++) {
            int[] pcsSel = model
                .component("comp1")
                .physics("ec")
                .feature("pcs" + pnum)
                .selection()
                .entities();
            if (pcsSel.length != 1) {
                System.out.println("Incorrect number of point selections for current");
                model.save(outpath + "/error_model.mph");
                System.exit(1);
            }
        }

        //set up study steps
        System.out.println("Running physics solutions");
        solutionSetup(model);

        for (int pnum = 1; pnum < numpcs + 1; pnum++) {
            //set all pcs to zero
            model.component("comp1").physics("ec").feature("pcs" + pnum).set("Qjp", "0");
        }
        model.save(outpath + "/debug_model.mph");

        //run solution
        for (int pnum = 0; pnum < numpcs; pnum++) {
            //set all pcs to zero
            model.component("comp1").physics("ec").feature("pcs" + (pnum + 1)).set("Qjp", ".001");
            model.sol("sol1").runAll(); //TODO: should really save these to their own directory (bases)
            System.out.println("Solution " + pnum + " solved, saving...");
            model.save(outpath + "/" + pnum + ".mph");
            System.out.println("Saved");
            model.component("comp1").physics("ec").feature("pcs" + (pnum + 1)).set("Qjp", "0");
        }

        System.out.println("All solutions solved successfully.");
        return model;
    }

    private static Model runDeformation(JSONObject config3D) throws IOException {
        Model model;
        // TODO: add error if cuff is not LivaNova
        ModelUtil.showProgress(true);
        // START PRE MESH
        System.out.println("Running Deformation Simulation...");
        // Define model object
        model = ModelUtil.createUnique("Model");

        // Define ModelWrapper class instance for model and projectPath
        String defdir = "doubledeform";
        String deformationConfigPath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("slides") +
            "/" +
            defdir +
            "/sipdefconfig.json";
        JSONObject defconfig = null;
        try {
            defconfig = JSONio.read(deformationConfigPath);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        model.param().set("para", "0", "ramping parameters");
        model.param().set("final_ratio", 0.5);
        model.param().set("id_start", "id_final1/final_ratio");
        assert defconfig != null;
        model
            .param()
            .set(
                "id_final1",
                (
                    defconfig.getJSONArray("epi_radius").getDouble(0) -
                    defconfig.getDouble("epipad")
                ) *
                2 +
                " [um]"
            );
        model
            .param()
            .set(
                "id_final2",
                (
                    defconfig.getJSONArray("epi_radius").getDouble(1) -
                    defconfig.getDouble("epipad")
                ) *
                2 +
                " [um]"
            );
        model.param().set("zcenter1", "5.2 [mm]"); //TODO: make this not hardcoded
        model.param().set("zcenter2", "zcenter1+8 [mm]"); //TODO: make this not hardcoded
        model.param().set("cufflen", "5.4 [mm]"); //TODO: make this not hardcoded
        model.param().set("sep", defconfig.getDouble("epipad") + "[um]");

        model.component().create("comp1", true);

        model.component("comp1").geom().create("geom1", 3);

        model.component().create("mcomp1", "MeshComponent");

        model.geom().create("mgeom1", 3);

        model.component("comp1").func().create("an1", "Analytic");
        model.component("comp1").func("an1").set("funcname", "rampFun");
        model.component("comp1").func("an1").set("expr", "(1-para)*2^(-para*10)");
        model.component("comp1").func("an1").set("args", new String[] { "para" });
        model.component("comp1").func("an1").set("argunit", new String[] { "1" });
        model
            .component("comp1")
            .func("an1")
            .set("plotargs", new String[][] { { "para", "0", "1" } });

        model.component("comp1").mesh().create("mesh3");
        model.mesh().create("mpart1", "mgeom1");

        model.geom("mgeom1").lengthUnit("mm");

        model.mesh("mpart1").create("imp1", "Import");
        model.mesh("mpart1").feature("imp1").set("source", "stl");
        model
            .mesh("mpart1")
            .feature("imp1")
            .set("filename", defconfig.getString("out") + "/predeform.stl");
        model.mesh("mpart1").feature("imp1").set("createdom", true);
        model.mesh("mpart1").feature("imp1").importData();
        model.mesh("mpart1").run();

        model.component("comp1").geom("geom1").lengthUnit("mm");
        model.component("comp1").geom("geom1").create("imp1", "Import");
        model.component("comp1").geom("geom1").feature("imp1").set("type", "mesh");
        model.component("comp1").geom("geom1").feature("imp1").set("mesh", "mpart1");
        model.component("comp1").geom("geom1").feature("imp1").set("defectremoval", 0.1);
        model.component("comp1").geom("geom1").feature("imp1").importData();
        model.component("comp1").geom("geom1").create("cyl1", "Cylinder");
        model
            .component("comp1")
            .geom("geom1")
            .feature("cyl1")
            .set("pos", new String[] { "0", "0", "zcenter1-cufflen/2" });
        model.component("comp1").geom("geom1").feature("cyl1").set("r", "id_start/2+0.5");
        model.component("comp1").geom("geom1").feature("cyl1").set("h", "cufflen");
        model.component("comp1").geom("geom1").create("cyl2", "Cylinder");
        model
            .component("comp1")
            .geom("geom1")
            .feature("cyl2")
            .set("pos", new String[] { "0", "0", "zcenter1-cufflen/2" });
        model.component("comp1").geom("geom1").feature("cyl2").set("r", "id_start/2");
        model.component("comp1").geom("geom1").feature("cyl2").set("h", "cufflen");
        model.component("comp1").geom("geom1").create("dif1", "Difference");
        model.component("comp1").geom("geom1").feature("dif1").selection("input").set("cyl1");
        model.component("comp1").geom("geom1").feature("dif1").selection("input2").set("cyl2");
        model.component("comp1").geom("geom1").create("tor1", "Torus");
        model
            .component("comp1")
            .geom("geom1")
            .feature("tor1")
            .set("pos", new String[] { "0", "0", "zcenter1+cufflen/2" });
        model.component("comp1").geom("geom1").feature("tor1").set("rmaj", "(id_start+0.5)/2");
        model.component("comp1").geom("geom1").feature("tor1").set("rmin", 0.25);
        model.component("comp1").geom("geom1").create("tor2", "Torus");
        model
            .component("comp1")
            .geom("geom1")
            .feature("tor2")
            .set("pos", new String[] { "0", "0", "zcenter1-cufflen/2" });
        model.component("comp1").geom("geom1").feature("tor2").set("rmaj", "(id_start+0.5)/2");
        model.component("comp1").geom("geom1").feature("tor2").set("rmin", 0.25);
        model.component("comp1").geom("geom1").create("uni1", "Union");
        model.component("comp1").geom("geom1").feature("uni1").set("intbnd", false);
        model
            .component("comp1")
            .geom("geom1")
            .feature("uni1")
            .selection("input")
            .set("dif1", "tor1", "tor2");
        model.component("comp1").geom("geom1").create("copy1", "Copy");
        model.component("comp1").geom("geom1").feature("copy1").set("displz", "8 [mm]");
        model.component("comp1").geom("geom1").feature("copy1").selection("input").set("uni1");
        model.component("comp1").geom("geom1").create("sca1", "Scale");
        model.component("comp1").geom("geom1").feature("sca1").set("type", "anisotropic");
        model
            .component("comp1")
            .geom("geom1")
            .feature("sca1")
            .set("factor", new String[] { "id_final2/id_final1", "id_final2/id_final1", "1" });
        model.component("comp1").geom("geom1").feature("sca1").selection("input").set("copy1");
        model.component("comp1").geom("geom1").feature("fin").label("Form Assembly");
        model.component("comp1").geom("geom1").feature("fin").set("action", "assembly");
        model.component("comp1").selection().create("box1", "Box");
        model.component("comp1").selection("box1").set("entitydim", 2);
        model.component("comp1").selection().create("adj1", "Adjacent");
        model.component("comp1").selection("adj1").set("entitydim", 2);
        model.component("comp1").selection("box1").set("inputent", "selections");
        model
            .component("comp1")
            .selection("box1")
            .set("input", new String[] { "geom1_imp1_mpart1_imp1_unnamed" });
        model.component("comp1").selection("box1").set("zmin", ".05");
        model.component("comp1").selection("box1").set("zmax", 18.35); //TODO make this not hardcoded
        model.component("comp1").selection("adj1").set("input", new String[] { "box1" });
        model.component("comp1").selection().create("com1", "Complement");
        model
            .component("comp1")
            .selection("com1")
            .set("input", new String[] { "geom1_imp1_mpart1_imp1_unnamed_1" });
        model.component("comp1").geom("geom1").run("fin");
        int[] out = model.component("comp1").selection("box1").entities();
        //raise error if length of out array is zero
        if (out.length == 0) {
            throw new RuntimeException(
                "No entities found in box selection. " +
                "Ensure that the geometry is correctly placed (bottom should be flush with z=0)."
            );
        }
        model.component("comp1").geom("geom1").create("cmf1", "CompositeFaces");
        model.component("comp1").geom("geom1").feature("cmf1").selection("input").set("fin", out);
        model.component("comp1").geom("geom1").run();

        model.component("comp1").pair().create("p1", "Contact");
        model
            .component("comp1")
            .pair("p1")
            .source()
            .set(
                7,
                8,
                9,
                10,
                11,
                12,
                16,
                17,
                18,
                19,
                20,
                21,
                31,
                32,
                33,
                34,
                35,
                36,
                40,
                41,
                42,
                43,
                44,
                45
            );
        model.component("comp1").pair("p1").destination().named("box1");
        model.component("comp1").pair().create("p2", "Contact");
        model.component("comp1").pair("p2").source().named("box1");
        model.component("comp1").pair("p2").destination().named("box1");

        model.component("comp1").material().create("mat1", "Common");
        model.component("comp1").material().create("mat2", "Common");
        model.component("comp1").material("mat1").selection().set(1, 2);
        model
            .component("comp1")
            .material("mat1")
            .propertyGroup()
            .create("Enu", "Young's modulus and Poisson's ratio");
        model
            .component("comp1")
            .material("mat2")
            .selection()
            .named("geom1_imp1_mpart1_imp1_unnamed_1");
        model
            .component("comp1")
            .material("mat2")
            .propertyGroup()
            .create("Enu", "Young's modulus and Poisson's ratio");

        model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");
        model.component("comp1").physics("solid").create("lemm2", "LinearElasticModel", 3);
        model.component("comp1").physics("solid").feature("lemm2").selection().set(1, 2);
        model
            .component("comp1")
            .physics("solid")
            .feature("lemm2")
            .create("eiel1", "ExternalStrain", 3);
        model.component("comp1").physics("solid").create("disp1", "Displacement3", 3);
        model.component("comp1").physics("solid").feature("disp1").selection().all();
        model.component("comp1").physics("solid").create("disp2", "Displacement2", 2);
        model.component("comp1").physics("solid").feature("disp2").selection().named("adj1");
        model.component("comp1").physics("solid").create("spf1", "SpringFoundation3", 3);
        model
            .component("comp1")
            .physics("solid")
            .feature("spf1")
            .selection()
            .named("geom1_imp1_mpart1_imp1_unnamed_1");
        model.component("comp1").physics("solid").create("disp3", "Displacement1", 1);
        model.component("comp1").physics("solid").feature("disp3").selection().set(19, 40, 67, 88);
        model.component("comp1").physics("solid").create("disp4", "Displacement1", 1);
        model.component("comp1").physics("solid").feature("disp4").selection().set(1, 48, 49, 96);
        model.component("comp1").physics("solid").create("cnt1", "SolidContact", 2);

        model.component("comp1").mesh("mesh3").create("ftri1", "FreeTri");
        model.component("comp1").mesh("mesh3").create("map1", "Map");
        model.component("comp1").mesh("mesh3").create("ftet1", "FreeTet");
        model.component("comp1").mesh("mesh3").feature("ftri1").create("size1", "Size");
        model.component("comp1").mesh("mesh3").feature("ftet1").create("size1", "Size");
        model.component("comp1").selection().create("sel1", "Explicit");
        model.component("comp1").selection("sel1").set("angletol", 5);
        model.component("comp1").selection("sel1").geom("geom1", 3, 2, new String[] { "exterior" });
        model.component("comp1").selection("sel1").set(1, 2);
        model.component("comp1").mesh("mesh3").feature("map1").selection().named("sel1");
        model.component("comp1").mesh("mesh3").feature("map1").create("size1", "Size");
        model.component("comp1").mesh("mesh3").feature("map1").feature("size1").set("hauto", 3);
        model.component("comp1").mesh("mesh3").feature("ftri1").selection().named("box1");
        model.component("comp1").pair("p1").searchMethod("direct");
        model.component("comp1").pair("p2").searchMethod("direct");
        model.component("comp1").pair("p2").manualDist(true);

        model.component("comp1").material("mat1").label("Alumina");
        model.component("comp1").material("mat1").set("family", "aluminum");
        model
            .component("comp1")
            .material("mat1")
            .propertyGroup("def")
            .set(
                "thermalexpansioncoefficient",
                new String[] { "8e-6[1/K]", "0", "0", "0", "8e-6[1/K]", "0", "0", "0", "8e-6[1/K]" }
            );
        model
            .component("comp1")
            .material("mat1")
            .propertyGroup("def")
            .set("heatcapacity", "900[J/(kg*K)]");
        model
            .component("comp1")
            .material("mat1")
            .propertyGroup("def")
            .set("density", "3900[kg/m^3]");
        model
            .component("comp1")
            .material("mat1")
            .propertyGroup("def")
            .set(
                "thermalconductivity",
                new String[] {
                    "27[W/(m*K)]",
                    "0",
                    "0",
                    "0",
                    "27[W/(m*K)]",
                    "0",
                    "0",
                    "0",
                    "27[W/(m*K)]",
                }
            );
        model.component("comp1").material("mat1").propertyGroup("Enu").set("E", "300e9[Pa]");
        model.component("comp1").material("mat1").propertyGroup("Enu").set("nu", "0.222");
        model.component("comp1").material("mat2").label("bendy");
        model.component("comp1").material("mat2").set("family", "aluminum");
        model
            .component("comp1")
            .material("mat2")
            .propertyGroup("def")
            .set(
                "thermalexpansioncoefficient",
                new String[] { "8e-6[1/K]", "0", "0", "0", "8e-6[1/K]", "0", "0", "0", "8e-6[1/K]" }
            );
        model
            .component("comp1")
            .material("mat2")
            .propertyGroup("def")
            .set("heatcapacity", "900[J/(kg*K)]");
        model
            .component("comp1")
            .material("mat2")
            .propertyGroup("def")
            .set("density", "3000[kg/m^3]");
        model
            .component("comp1")
            .material("mat2")
            .propertyGroup("def")
            .set(
                "thermalconductivity",
                new String[] {
                    "27[W/(m*K)]",
                    "0",
                    "0",
                    "0",
                    "27[W/(m*K)]",
                    "0",
                    "0",
                    "0",
                    "27[W/(m*K)]",
                }
            );
        model.component("comp1").material("mat2").propertyGroup("Enu").set("E", "100e3[Pa]");
        model.component("comp1").material("mat2").propertyGroup("Enu").set("nu", "0.3");
        model
            .component("comp1")
            .physics("solid")
            .prop("AdvancedSettings")
            .set("GroupPhysOdesRc", false);
        model
            .component("comp1")
            .physics("solid")
            .prop("AdvancedSettings")
            .set("GroupPhysOdesAtt", false);
        model
            .component("comp1")
            .physics("solid")
            .feature("dcnt1")
            .set("ContactMethodCtrl", "Nitsche");
        model.component("comp1").physics("solid").feature("lemm2").label("compressor");
        model
            .component("comp1")
            .physics("solid")
            .feature("lemm2")
            .feature("eiel1")
            .set(
                "stch",
                new String[][] {
                    { "1 - (1-final_ratio)*(1-rampFun(para))" },
                    { "1 - (1-final_ratio)*(1-rampFun(para))" },
                    { "1" },
                }
            );
        model
            .component("comp1")
            .physics("solid")
            .feature("lemm2")
            .feature("eiel1")
            .set("StrainInput", "Stretches");
        model
            .component("comp1")
            .physics("solid")
            .feature("disp1")
            .set("Direction", new int[][] { { 0 }, { 0 }, { 1 } });
        model.component("comp1").physics("solid").feature("disp1").label("fix z");
        model
            .component("comp1")
            .physics("solid")
            .feature("disp2")
            .set("Direction", new int[][] { { 1 }, { 1 }, { 1 } });
        model.component("comp1").physics("solid").feature("disp2").label("fix ends");
        model
            .component("comp1")
            .physics("solid")
            .feature("spf1")
            .set(
                "kPerVolume",
                new String[][] {
                    { "1e4*rampFun(para)" },
                    { "0" },
                    { "0" },
                    { "0" },
                    { "1e4*rampFun(para)" },
                    { "0" },
                    { "0" },
                    { "0" },
                    { "1e4*rampFun(para)" },
                }
            );
        model
            .component("comp1")
            .physics("solid")
            .feature("disp3")
            .set("Direction", new int[][] { { 1 }, { 0 }, { 0 } });
        model.component("comp1").physics("solid").feature("disp3").label("fix x");
        model
            .component("comp1")
            .physics("solid")
            .feature("disp4")
            .set("Direction", new int[][] { { 0 }, { 1 }, { 0 } });
        model.component("comp1").physics("solid").feature("disp4").label("fix y");
        model
            .component("comp1")
            .physics("solid")
            .feature("cnt1")
            .set("ContactMethodCtrl", "Nitsche");
        model.component("comp1").physics("solid").feature("cnt1").set("pairs", "p2");

        model.component("comp1").mesh("mesh3").feature("size").set("hauto", 4);
        model.component("comp1").mesh("mesh3").feature("ftri1").feature("size1").set("hauto", 3);
        model.component("comp1").mesh("mesh3").feature("ftet1").feature("size1").set("hauto", 5);
        model.save(defconfig.getString("out") + "/premesh_def_debug.mph");
        model.component("comp1").mesh("mesh3").run();

        model.study().create("std1");
        model.study("std1").create("stat", "Stationary");

        model.sol().create("sol1");
        model.sol("sol1").study("std1");
        model.sol("sol1").attach("std1");
        model.sol("sol1").create("st1", "StudyStep");
        model.sol("sol1").create("v1", "Variables");
        model.sol("sol1").create("s1", "Stationary");
        model.sol("sol1").feature("s1").create("p1", "Parametric");
        model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
        model.sol("sol1").feature("s1").create("d1", "Direct");
        model.sol("sol1").feature("s1").create("i1", "Iterative");
        model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("pr")
            .create("so1", "SOR");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("po")
            .create("so1", "SOR");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .create("mg1", "Multigrid");
        model.sol("sol1").feature("s1").feature().remove("fcDef");

        model.study("std1").feature("stat").set("plot", true);
        model.study("std1").feature("stat").set("useparam", true);
        model.study("std1").feature("stat").set("pname", new String[] { "para" });
        model
            .study("std1")
            .feature("stat")
            .set("plistarr", new String[] { "1e-3 range(1e-2 ,1e-2,1)" });
        model.study("std1").feature("stat").set("punit", new String[] { "" });

        model.sol("sol1").attach("std1");
        model.sol("sol1").feature("st1").label("Compile Equations: Stationary");
        model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
        model.sol("sol1").feature("v1").set("clistctrl", new String[] { "p1" });
        model.sol("sol1").feature("v1").set("cname", new String[] { "para" });
        model.sol("sol1").feature("v1").set("clist", new String[] { "1e-3 range(1e-2 ,1e-2,1)" });
        model.sol("sol1").feature("v1").feature("comp1_u").set("scalemethod", "manual");
        model
            .sol("sol1")
            .feature("v1")
            .feature("comp1_u")
            .set("scaleval", "1e-2*0.011204025887081638");
        model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
        model.sol("sol1").feature("s1").set("probesel", "none");
        model.sol("sol1").feature("s1").feature("dDef").label("Direct 2");
        model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
        model.sol("sol1").feature("s1").feature("aDef").set("cachepattern", true);
        model.sol("sol1").feature("s1").feature("p1").label("Parametric 1.1");
        model.sol("sol1").feature("s1").feature("p1").set("pname", new String[] { "para" });
        model
            .sol("sol1")
            .feature("s1")
            .feature("p1")
            .set("plistarr", new String[] { "1e-3 range(1e-2 ,1e-2,1)" });
        model.sol("sol1").feature("s1").feature("p1").set("punit", new String[] { "" });
        model.sol("sol1").feature("s1").feature("p1").set("porder", "constant");
        model.sol("sol1").feature("s1").feature("p1").set("uselsqdata", false);
        model.sol("sol1").feature("s1").feature("p1").set("plot", true);
        model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
        model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "d1");
        model.sol("sol1").feature("s1").feature("fc1").set("dtech", "ddog");
        model.sol("sol1").feature("s1").feature("d1").label("Suggested Direct Solver (solid)");
        model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
        model.sol("sol1").feature("s1").feature("d1").set("pivotperturb", 1.0E-9);
        model.sol("sol1").feature("s1").feature("d1").set("nliniterrefine", true);
        model.sol("sol1").feature("s1").feature("i1").label("Suggested Iterative Solver (solid)");
        model.sol("sol1").feature("s1").feature("i1").set("prefuntype", "right");
        model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", true);
        model.sol("sol1").feature("s1").feature("i1").set("rhob", 40);
        model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
        model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("pr")
            .label("Presmoother 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("pr")
            .feature("soDef")
            .label("SOR 2");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("pr")
            .feature("so1")
            .label("SOR 1.1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("pr")
            .feature("so1")
            .set("iter", 1);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("po")
            .label("Postsmoother 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("po")
            .feature("soDef")
            .label("SOR 2");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("po")
            .feature("so1")
            .label("SOR 1.1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("po")
            .feature("so1")
            .set("iter", 1);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .label("Coarse Solver 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("dDef")
            .label("Direct 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .label("Multigrid 1.1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .set("prefun", "saamg");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .set("iter", 2);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .set("mglevels", 2);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .set("maxcoarsedof", 10000);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .set("usesmooth", false);
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("pr")
            .label("Presmoother 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("pr")
            .feature("soDef")
            .label("SOR 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("po")
            .label("Postsmoother 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("po")
            .feature("soDef")
            .label("SOR 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("cs")
            .label("Coarse Solver 1");
        model
            .sol("sol1")
            .feature("s1")
            .feature("i1")
            .feature("mg1")
            .feature("cs")
            .feature("mg1")
            .feature("cs")
            .feature("dDef")
            .label("Direct 1");
        model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 5);
        model.save(defconfig.getString("out") + "/deformation_debug.mph");
        System.out.println(
            "temporary exit, manually solve deformation. Remember to look over the sip input."
        );
        System.exit(0);
        model.sol("sol1").runAll();
        model.result().dataset("dset1").createDeformedConfig("deform1", "mesh4");
        model.component("comp1").geom("deform1").export().setType("stlbin");
        model
            .component("comp1")
            .geom("deform1")
            .export(defconfig.getString("out") + "/postdeform.stl");

        model.result().create("pg1", "PlotGroup3D");
        model.result().create("pg2", "PlotGroup3D");
        model.result().create("pg3", "PlotGroup3D");
        model.result("pg1").create("vol1", "Volume");
        model.result("pg1").feature("vol1").set("expr", "solid.mises");
        model.result("pg1").feature("vol1").create("def", "Deform");
        model.result("pg2").create("vol1", "Volume");
        model.result("pg2").feature("vol1").create("def", "Deform");
        model.result("pg3").create("vol1", "Volume");
        model.result("pg3").feature("vol1").set("expr", "solid.mises");
        model.result("pg3").feature("vol1").create("def", "Deform");
        model.result("pg1").label("Stress (solid)");
        model.result("pg1").set("frametype", "spatial");
        model
            .result("pg1")
            .feature("vol1")
            .set(
                "const",
                new String[][] {
                    {
                        "solid.refpntx",
                        "0",
                        "Reference point for moment computation, x-coordinate",
                    },
                    {
                        "solid.refpnty",
                        "0",
                        "Reference point for moment computation, y-coordinate",
                    },
                    {
                        "solid.refpntz",
                        "0",
                        "Reference point for moment computation, z-coordinate",
                    },
                }
            );
        model.result("pg1").feature("vol1").set("colortable", "Prism");
        model.result("pg1").feature("vol1").set("colorscalemode", "logarithmic");
        model.result("pg1").feature("vol1").set("resolution", "custom");
        model.result("pg1").feature("vol1").set("refine", 2);
        model.result("pg1").feature("vol1").set("threshold", "manual");
        model.result("pg1").feature("vol1").set("thresholdvalue", 0.2);
        model.result("pg1").feature("vol1").set("resolution", "custom");
        model.result("pg1").feature("vol1").set("refine", 2);
        model.result("pg1").feature("vol1").feature("def").set("scaleactive", true);
        model.result("pg2").label("Stress (solid) 1");
        model.result("pg2").set("looplevel", new int[] { 52 });
        model.result("pg2").set("frametype", "spatial");
        model
            .result("pg2")
            .feature("vol1")
            .set(
                "const",
                new String[][] {
                    {
                        "solid.refpntx",
                        "0",
                        "Reference point for moment computation, x-coordinate",
                    },
                    {
                        "solid.refpnty",
                        "0",
                        "Reference point for moment computation, y-coordinate",
                    },
                    {
                        "solid.refpntz",
                        "0",
                        "Reference point for moment computation, z-coordinate",
                    },
                }
            );
        model.result("pg2").feature("vol1").set("colortable", "Prism");
        model.result("pg2").feature("vol1").set("resolution", "custom");
        model.result("pg2").feature("vol1").set("refine", 2);
        model.result("pg2").feature("vol1").set("threshold", "manual");
        model.result("pg2").feature("vol1").set("thresholdvalue", 0.2);
        model.result("pg2").feature("vol1").set("resolution", "custom");
        model.result("pg2").feature("vol1").set("refine", 2);
        model.result("pg2").feature("vol1").feature("def").set("scaleactive", true);
        model.result("pg3").label("Stress (solid) 2");
        model.result("pg3").set("frametype", "spatial");
        model
            .result("pg3")
            .feature("vol1")
            .set(
                "const",
                new String[][] {
                    {
                        "solid.refpntx",
                        "0",
                        "Reference point for moment computation, x-coordinate",
                    },
                    {
                        "solid.refpnty",
                        "0",
                        "Reference point for moment computation, y-coordinate",
                    },
                    {
                        "solid.refpntz",
                        "0",
                        "Reference point for moment computation, z-coordinate",
                    },
                }
            );
        model.result("pg3").feature("vol1").set("colortable", "Prism");
        model.result("pg3").feature("vol1").set("colorscalemode", "logarithmic");
        model.result("pg3").feature("vol1").set("resolution", "custom");
        model.result("pg3").feature("vol1").set("refine", 2);
        model.result("pg3").feature("vol1").set("threshold", "manual");
        model.result("pg3").feature("vol1").set("thresholdvalue", 0.2);
        model.result("pg3").feature("vol1").set("resolution", "custom");
        model.result("pg3").feature("vol1").set("refine", 2);
        model.result("pg3").feature("vol1").feature("def").set("scaleactive", true);
        model.result().dataset("dset1").selection().geom("geom1", 3);
        model.result().dataset("dset1").selection().named("geom1_imp1_mpart1_imp1_unnamed_1");
        model.result("pg1").run();
        model.save(defconfig.getString("out") + "/deformation.mph");
        return model;
    }

    private static Model buildGeometry(
        JSONObject config3D,
        JSONObject modelData,
        String projectPath
    ) throws FileNotFoundException {
        ModelWrapper mw;
        Model model;
        String mediumPrimitiveString = "Medium_Primitive";
        String instanceLabelDistalMedium = DISTAL_MEDIUM;
        String instanceLabelProximalMedium = PROXIMAL_MEDIUM;
        String outpath =
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("comsol");

        String geomFile = outpath + "/geometry.mph";

        // START PRE MESH
        System.out.println("Generating Geometry...");
        // Define model object
        model = ModelUtil.createUnique("Model");
        // Add component node 1
        model.component().create("comp1", true);
        // Add 3D geom to component node 1
        model.component("comp1").geom().create("geom1", 3);
        // geometry shape order
        assert modelData != null;
        String sorder = modelData.getJSONObject("solver").getString("sorder");
        model.component("comp1").sorder(sorder);
        // Set default length units to micron
        model.component("comp1").geom("geom1").lengthUnit("\u00b5m");
        // Add materials node to component node 1
        model.component("comp1").physics().create("ec", "ConductiveMedia", "geom1");
        // and mesh node to component node 1
        model.component("comp1").mesh().create("mesh1");

        // Define ModelWrapper class instance for model and projectPath
        mw = new ModelWrapper(model, projectPath);

        // FEM MODEL GEOMETRY
        // Set MEDIUM parameters
        JSONObject distalMedium = modelData.getJSONObject("medium").getJSONObject("distal");
        JSONObject proximalMedium = modelData.getJSONObject("medium").getJSONObject("proximal");

        String mediumParamsLabel = "Medium Parameters";
        ModelParamGroup mediumParams = model.param().group().create(mediumParamsLabel);
        mediumParams.label(mediumParamsLabel);

        double proximal_length = proximalMedium.getDouble("length");
        double proximal_radius = proximalMedium.getDouble("radius");

        String bounds_unit = "[um]";
        mediumParams.set("z_nerve", proximal_length + " " + bounds_unit);
        mediumParams.set("r_proximal", proximal_radius + " " + bounds_unit);

        double prox_z = proximalMedium.getJSONObject("shift").getDouble("z");

        mediumParams.set("prox_shift_x", 0 + " " + bounds_unit);
        mediumParams.set("prox_shift_y", 0 + " " + bounds_unit);
        mediumParams.set("prox_shift_z", prox_z + " " + bounds_unit);

        if (distalMedium.getBoolean("exist")) {
            double distal_length = distalMedium.getDouble("length");
            double distal_radius = distalMedium.getDouble("radius");
            double distal_x = 0;
            double distal_y = 0;
            double distal_z = distalMedium.getJSONObject("shift").getDouble("z");

            mediumParams.set("z_distal", distal_length + " " + bounds_unit);
            mediumParams.set("r_distal", distal_radius + " " + bounds_unit);
            mediumParams.set("distal_shift_x", distal_x + " " + bounds_unit);
            mediumParams.set("distal_shift_y", distal_y + " " + bounds_unit);
            mediumParams.set("distal_shift_z", distal_z + " " + bounds_unit);
        }

        // Create PART PRIMITIVE for MEDIUM
        String partID = mw.im.next("part", mediumPrimitiveString);
        IdentifierManager partPrimitiveIM;
        try {
            partPrimitiveIM =
                Part.createEnvironmentPartPrimitive(partID, mediumPrimitiveString, mw);
            mw.partPrimitiveIMs.put(mediumPrimitiveString, partPrimitiveIM);
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        }

        // Create PART INSTANCES for MEDIUM (Distal and Proximal)
        if (distalMedium.getBoolean("exist")) {
            String mediumDistal_instanceID = mw.im.next("pi", instanceLabelDistalMedium);

            if (proximalMedium.getBoolean("distant_ground")) {
                System.out.println(
                    "WARNING: you have a distal domain, as well as a proximal domain " +
                    "that is grounded... make sure this is something you actually want to do..."
                );
            }

            try {
                Part.createEnvironmentPartInstance(
                    mediumDistal_instanceID,
                    instanceLabelDistalMedium,
                    mediumPrimitiveString,
                    mw,
                    distalMedium
                );
            } catch (IllegalArgumentException e) {
                e.printStackTrace();
            }
        }

        String mediumProximal_instanceID = mw.im.next("pi", instanceLabelProximalMedium);
        try {
            Part.createEnvironmentPartInstance(
                mediumProximal_instanceID,
                instanceLabelProximalMedium,
                mediumPrimitiveString,
                mw,
                proximalMedium
            );
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        }

        ModelParamGroup nerveParams;
        // Set NERVE MORPHOLOGY parameters
        String morphology_unit = "um";
        String nerveParamsLabal = "Nerve Parameters";
        nerveParams = model.param().group().create(nerveParamsLabal);
        nerveParams.label(nerveParamsLabal);

        nerveParams.set(
            "nerve_length",
            modelData.getDouble("nerve_length") + " [" + morphology_unit + "]"
        );
        nerveParams.set(
            "r_nerve1",
            modelData.getDouble("min_radius_enclosing_circle1") + " [" + morphology_unit + "]"
        );
        nerveParams.set(
            "r_nerve2",
            modelData.getDouble("min_radius_enclosing_circle2") + " [" + morphology_unit + "]"
        );

        // Set CUFF POSITIONING parameters
        String cuffConformationParamsLabel = "Cuff Conformation Parameters";
        ModelParamGroup cuffConformationParams = model
            .param()
            .group()
            .create(cuffConformationParamsLabel);
        cuffConformationParams.label(cuffConformationParamsLabel);

        String cuff_shift_unit = "[micrometer]";
        String cuff_rot_unit = "[degree]";
        Double cuff_shift_x = modelData.getJSONObject("cuff").getJSONObject("shift").getDouble("x");
        Double cuff_shift_y = modelData.getJSONObject("cuff").getJSONObject("shift").getDouble("y");
        Double cuff_shift_z = modelData.getJSONObject("cuff").getJSONObject("shift").getDouble("z");
        Double cuff_rot_pos = modelData
            .getJSONObject("cuff")
            .getJSONObject("rotate")
            .getDouble("pos_ang");
        Double cuff_rot_add = modelData
            .getJSONObject("cuff")
            .getJSONObject("rotate")
            .getDouble("add_ang");

        cuffConformationParams.set("cuff_shift_x", cuff_shift_x + " " + cuff_shift_unit);
        cuffConformationParams.set("cuff_shift_y", cuff_shift_y + " " + cuff_shift_unit);
        cuffConformationParams.set("cuff_shift_z", cuff_shift_z + " " + cuff_shift_unit);
        cuffConformationParams.set("cuff_rot", cuff_rot_pos + cuff_rot_add + " " + cuff_rot_unit);

        // add PART PRIMITIVES for CUFF
        // Read cuff to build from model.json (cuff.preset) which links to JSON containing instantiations of parts
        JSONObject cuffObject = (JSONObject) modelData.get("cuff");
        String cuff = cuffObject.getString("preset");
        mw.addCuffPartPrimitives(cuff);

        try {
            System.out.println("Saving pre-run geometry MPH file to: " + geomFile);
            model.save(geomFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        // add PART INSTANCES for cuff
        mw.addCuffPartInstances(
            cuff,
            config3D.getString("project_path") +
            "/" +
            config3D.getJSONObject("path").getString("comsol") +
            "/pcs"
        );

        try {
            System.out.println("Saving pre-run geometry MPH file to: " + geomFile);
            model.save(geomFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // create UNIONS
        mw.createUnions();

        model.component("comp1").geom("geom1").run("fin");

        //TODO build mesh

        String stlpath = outpath + "/stl";
        File stlPathFile = new File(stlpath);
        if (!stlPathFile.exists()) {
            boolean success = stlPathFile.mkdirs();
            assert success;
        }

        System.out.println("Exporting geometry STL files to: " + stlpath);

        //export STLs
        JSONObject cuffData = JSONio.read(
            String.join("/", new String[] { mw.root, "config", "system", "cuffs", cuff })
        );

        HashMap<String, String> matmap = new HashMap<>();
        JSONObject dommap = new JSONObject();
        //create dummy material for doing selections
        model.component("comp1").material().create("mat1", "Common");
        //run mesh
        model.component("comp1").mesh("mesh1").autoMeshSize(3); //set mesh finer for more complex geometries
        model.component("comp1").mesh("mesh1").run();
        //export as nastran
        model.component("comp1").mesh("mesh1").export().set("type", "nastran");
        model.component("comp1").mesh("mesh1").export(stlpath + "/alldomain.nas");

        // cuff stl export
        for (Object item : (JSONArray) cuffData.get("instances")) {
            JSONObject itemObject = (JSONObject) item;

            String instanceLabel = (String) itemObject.get("label");
            String pseudonym = (String) itemObject.get("type");

            IdentifierManager myIM = mw.getPartPrimitiveIM(pseudonym);
            if (myIM == null) throw new IllegalArgumentException(
                "IdentfierManager not created for name: " + pseudonym
            );

            String[] myLabels = myIM.labels; // may be null, but that is ok if not used

            // assign cuff materials
            JSONArray materials = itemObject.getJSONArray("materials");
            for (Object o : materials) {
                int label_index = ((JSONObject) o).getInt("label_index");
                String selection = myLabels[label_index];
                String info = ((JSONObject) o).getString("info");
                if (info.equals("fill")) { //maybe skip this
                    model.component("comp1").geom("geom1").export().selection().init();
                    String thislabel = mw.im.get(instanceLabel);
                    model.component("comp1").geom("geom1").export().selection().set(thislabel);
                } else {
                    model.component("comp1").geom("geom1").export().selection().init(3);
                    String thislabel = mw.im.get(instanceLabel) + "_" + myIM.get(selection);
                    model.component("comp1").geom("geom1").export().selection().named(thislabel);
                }
                if (
                    info.equals("recess") &
                    !config3D.getJSONObject("mesh").getBoolean("use_nastran")
                ) {
                    System.out.println("Skipping recess export, not implemented");
                } else {
                    model.component("comp1").geom("geom1").export().setSTLFormat("binary");
                    model
                        .component("comp1")
                        .geom("geom1")
                        .export(stlpath + "/" + instanceLabel + "-" + selection + ".stl");
                    matmap.put(instanceLabel + "-" + selection, info);
                }
                // Need to record the domains associated with each cuff material. First set mat1 selection
                //example: model.component("comp1").material("mat1").selection().named("geom1_pi6_csel7_dom");
                try {
                    model
                        .component("comp1")
                        .material("mat1")
                        .selection()
                        .named(
                            "geom1_" + mw.im.get(instanceLabel) + "_" + myIM.get(selection) + "_dom"
                        );
                    int[] ents = model.component("comp1").material("mat1").selection().entities(3);
                    dommap.put(instanceLabel + "-" + selection, ents);
                } catch (Exception e) {
                    System.out.println(
                        "No domain associated with this cuff material:" +
                        instanceLabel +
                        "-" +
                        selection
                    );
                }
            }
        }
        model.component("comp1").material().remove("mat1");
        // medium stl export
        model.component("comp1").geom("geom1").export().setSTLFormat("binary");
        model.component("comp1").geom("geom1").export().selection().init(2);
        model
            .component("comp1")
            .geom("geom1")
            .export()
            .selection()
            .named(mw.im.get(instanceLabelProximalMedium) + "_csel1");
        model.component("comp1").geom("geom1").export(stlpath + "/medium.stl");
        matmap.put("medium", "medium");

        JSONObject matmapout = new JSONObject(matmap);
        JSONio.write(stlpath + "/matmap.json", matmapout);
        JSONio.write(stlpath + "/dommap.json", dommap);

        // end export stl
        // save IM !!!!

        String imFile = outpath + "/im.json";
        JSONio.write(imFile, mw.im.toJSONObject()); // write to file

        // Saving model post-run geometry for debugging
        try {
            System.out.println("Saving geometry MPH file to: " + geomFile);
            model.save(geomFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return model;
    }

    private static void solutionSetup(Model model) {
        model.study().create("std1");
        model.study("std1").create("stat", "Stationary");
        model.study("std1").feature("stat").activate("ec", true);
        model.study("std1").feature("stat").activate("cc", true);
        model.study("std1").create("stat2", "Stationary");
        model.study("std1").feature("stat").setIndex("activate", false, 1);
        model.study("std1").feature("stat2").setIndex("activate", false, 3);

        model.sol().create("sol1");
        model.sol("sol1").study("std1");

        model.study("std1").feature("stat").set("notlistsolnum", 1);
        model.study("std1").feature("stat").set("notsolnum", "1");
        model.study("std1").feature("stat").set("listsolnum", 1);
        model.study("std1").feature("stat").set("solnum", "1");
        model.study("std1").feature("stat2").set("notlistsolnum", 1);
        model.study("std1").feature("stat2").set("notsolnum", "auto");
        model.study("std1").feature("stat2").set("listsolnum", 1);
        model.study("std1").feature("stat2").set("solnum", "auto");

        model.sol("sol1").create("st1", "StudyStep");
        model.sol("sol1").feature("st1").set("study", "std1");
        model.sol("sol1").feature("st1").set("studystep", "stat");
        model.sol("sol1").create("v1", "Variables");
        model.sol("sol1").feature("v1").set("control", "stat");
        model.sol("sol1").create("s1", "Stationary");
        model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
        model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "dDef");
        model.sol("sol1").feature("s1").feature().remove("fcDef");
        model.sol("sol1").create("su1", "StoreSolution");
        model.sol("sol1").create("st2", "StudyStep");
        model.sol("sol1").feature("st2").set("study", "std1");
        model.sol("sol1").feature("st2").set("studystep", "stat2");
        model.sol("sol1").create("v2", "Variables");
        model.sol("sol1").feature("v2").set("initmethod", "sol");
        model.sol("sol1").feature("v2").set("initsol", "sol1");
        model.sol("sol1").feature("v2").set("notsolmethod", "sol");
        model.sol("sol1").feature("v2").set("notsol", "sol1");
        model.sol("sol1").feature("v2").set("control", "stat2");
        model.sol("sol1").create("s2", "Stationary");
        model.sol("sol1").feature("s2").create("fc1", "FullyCoupled");
        model.sol("sol1").feature("s2").create("i1", "Iterative");
        model.sol("sol1").feature("s2").feature("i1").set("linsolver", "cg");
        model.sol("sol1").feature("s2").feature("i1").create("mg1", "Multigrid");
        model.sol("sol1").feature("s2").feature("i1").feature("mg1").set("prefun", "amg");
        model.sol("sol1").feature("s2").feature("fc1").set("linsolver", "i1");
        model.sol("sol1").feature("s2").feature().remove("fcDef");
        model.sol("sol1").feature("v2").set("notsolnum", "auto");
        model.sol("sol1").feature("v2").set("notsolvertype", "solnum");
        model.sol("sol1").feature("v2").set("notlistsolnum", new String[] { "1" });
        model.sol("sol1").feature("v2").set("notsolnum", "auto");
        model.sol("sol1").feature("v2").set("notlistsolnum", new String[] { "1" });
        model.sol("sol1").feature("v2").set("notsolnum", "auto");
        model.sol("sol1").feature("v2").set("control", "stat2");
        model.sol("sol1").feature("v2").set("solnum", "auto");
        model.sol("sol1").feature("v2").set("solvertype", "solnum");
        model.sol("sol1").feature("v2").set("listsolnum", new String[] { "1" });
        model.sol("sol1").feature("v2").set("solnum", "auto");
        model.sol("sol1").feature("v2").set("listsolnum", new String[] { "1" });
        model.sol("sol1").feature("v2").set("solnum", "auto");
        model.sol("sol1").attach("std1");
    }

    private static void LicenseCheckout(long wait_hours, String licenseType)
        throws InterruptedException {
        System.out.println(
            "Attempting to check out " +
            licenseType +
            " license. System will wait up to " +
            wait_hours +
            " hours for an available license seat."
        );
        boolean lic = false;
        long start = System.currentTimeMillis();
        long stop = wait_hours * 60 * 60 * 1000 + start;
        while (System.currentTimeMillis() < stop) {
            lic = ModelUtil.checkoutLicense(licenseType);
            if (lic) {
                long now = System.currentTimeMillis();
                double elapsed =
                    (Long.valueOf(now).doubleValue() - Long.valueOf(start).doubleValue()) /
                    (60 * 60 * 1000);
                System.out.printf(
                    licenseType + " license seat obtained (took %.3f hours).%n",
                    elapsed
                );
                break;
            } else {
                TimeUnit.SECONDS.sleep(600);
            }
        }
        if (!lic) {
            System.out.println(
                licenseType +
                " license did not become available within the specified time window. Exiting..."
            );
            System.exit(1);
        }
    }
}
