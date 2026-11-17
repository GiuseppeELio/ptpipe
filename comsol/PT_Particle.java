/*
 * PT_Particle.java  -  ptpipe COMSOL template, levels L1 + L3 (particle scale, nm)
 *
 * One or more Au nanoparticles in a homogeneous medium:
 *   shape = sphere  : sphere of radius R_np at the origin
 *   shape = rod     : prolate spheroid (semi-axes a_rod, b_rod, b_rod), long axis along x
 *   shape = cluster : spheres read from a text file (x y z R in nm, one particle per line)
 *
 * 1) Electromagnetic Waves, Frequency Domain, scattered-field formulation, spherical PML.
 *    sigma_abs = int(Qh)/I0 over the particles; sigma_sca = Poynting flux of the scattered
 *    field through an integration sphere (r_sca); sigma_ext = sum.
 *    pol = 1: k||z, E||x    pol = 2: k||z, E||y    pol = 3: k||x, E||z
 * 2) Heat Transfer in Solids in the EM domain + an outer thermal shell (T = T0 at R_th):
 *    source = emw.Qh * I_phys/I0 inside the particles (linear problem).
 *
 * Particle material (key "material"): "Au" = analytic model of Etchegoin et al., J. Chem. Phys.
 * 125, 164705 (2006), identical to ptpipe.materials.eps_au_etchegoin; "nk:materials/<file>" =
 * tabulated n,k (lambda_nm n k); or a real index. The same data feed the Python references, so
 * COMSOL vs Python compares solvers, not data.
 *
 * Configuration: comsol/particle.cfg (key = value; see comsol/configs/).
 * Output: results/<tag>_*.txt.   Run: comsol/run_validation.sh
 * ROOT: path of the repository; edit the constant below or pass -Dptpipe.root=<path>.
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class PT_Particle {

  // repository root (trailing slash). Edit once, or pass -Dptpipe.root=<path> to the JVM.
  static String BASE = "/path/to/ptpipe/";   // set automatically by run_validation.sh

  public static Model run() {
    try { return runImpl(); } catch (Exception e) { throw new RuntimeException(e); }
  }

  static Map<String, String> readCfg(String f) throws IOException {
    Map<String, String> c = new LinkedHashMap<>();
    BufferedReader br = new BufferedReader(new FileReader(f));
    String s;
    while ((s = br.readLine()) != null) {
      int h = s.indexOf('#'); if (h >= 0) s = s.substring(0, h);
      int e = s.indexOf('='); if (e < 0) continue;
      c.put(s.substring(0, e).trim(), s.substring(e + 1).trim());
    }
    br.close();
    return c;
  }
  static String g(Map<String, String> c, String k, String d) { return c.containsKey(k) ? c.get(k) : d; }
  static double[] list(String s) {
    s = s.trim();
    if (s.contains(":")) {             // a:step:b
      String[] p = s.split(":");
      double a = Double.parseDouble(p[0]), st = Double.parseDouble(p[1]), b = Double.parseDouble(p[2]);
      int n = (int) Math.floor((b - a) / st + 1e-9) + 1;
      double[] o = new double[n]; for (int i = 0; i < n; i++) o[i] = a + i * st; return o;
    }
    String[] p = s.split(",");
    double[] o = new double[p.length]; for (int i = 0; i < p.length; i++) o[i] = Double.parseDouble(p[i].trim()); return o;
  }

  static Model runImpl() throws Exception {
    try { String r = System.getProperty("ptpipe.root"); if (r != null) BASE = r.endsWith("/") ? r : r + "/"; } catch (Exception e) { }
    Map<String, String> cfg = readCfg(BASE + "comsol/particle.cfg");
    String tag = g(cfg, "tag", "particle");
    String shape = g(cfg, "shape", "sphere");
    PrintWriter log = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_log.txt"), true);
    log.println("cfg: " + cfg);

    double[] lams = list(g(cfg, "lams", "450:10:700"));
    double[] pols = list(g(cfg, "pols", "1"));
    double lamMin = Arrays.stream(lams).min().getAsDouble(), lamMax = Arrays.stream(lams).max().getAsDouble();
    double lamHeat = Double.parseDouble(g(cfg, "lam_heat", "530"));
    int polHeat = (int) Double.parseDouble(g(cfg, "pol_heat", "1"));
    lamMin = Math.min(lamMin, lamHeat); lamMax = Math.max(lamMax, lamHeat);

    // particles: x y z R (nm)
    List<double[]> P = new ArrayList<>();
    double R_np = Double.parseDouble(g(cfg, "R_np", "20"));
    double a_rod = Double.parseDouble(g(cfg, "a_rod", "40")), b_rod = Double.parseDouble(g(cfg, "b_rod", "10"));
    if (shape.equals("cluster")) {
      BufferedReader br = new BufferedReader(new FileReader(BASE + "comsol/" + g(cfg, "cluster_file", "cluster.txt")));
      String s;
      while ((s = br.readLine()) != null) {
        s = s.trim(); if (s.isEmpty() || s.startsWith("#") || s.startsWith("%")) continue;
        String[] t = s.split("[\\s,]+");
        P.add(new double[]{Double.parseDouble(t[0]), Double.parseDouble(t[1]), Double.parseDouble(t[2]), Double.parseDouble(t[3])});
      }
      br.close();
    } else if (shape.equals("rod")) {
      P.add(new double[]{0, 0, 0, a_rod});
    } else {
      P.add(new double[]{0, 0, 0, R_np});
    }
    int N = P.size();
    double Rext = 0, Rmin = 1e9;
    for (double[] p : P) { Rext = Math.max(Rext, Math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) + p[3]); Rmin = Math.min(Rmin, p[3]); }
    double hNp = shape.equals("rod") ? b_rod / 4 : Rmin / 4;
    if (cfg.containsKey("h_np")) hNp = Double.parseDouble(cfg.get("h_np"));
    log.printf("N = %d particles, R_ext = %.2f nm, h_np = %.2f nm%n", N, Rext, hNp);

    Model model = ModelUtil.create("Model");
    model.modelPath(BASE + "comsol");
    model.label("PT_Particle_" + tag + ".mph");

    // ------------------------------------------------------------------ parameters
    model.param().set("lam0", lamHeat + "[nm]", "Wavelength (swept in the spectrum)");
    model.param().set("lam_nm", "lam0/1[nm]");
    model.param().set("lam_min", lamMin + "[nm]");
    model.param().set("lam_max", lamMax + "[nm]");
    model.param().set("pol", String.valueOf(polHeat), "1: k||z E||x; 2: k||z E||y; 3: k||x E||z");
    model.param().set("n_m", g(cfg, "n_m", "1.33"), "Refractive index of the medium");
    model.param().set("R_np", R_np + "[nm]");
    model.param().set("a_rod", a_rod + "[nm]");
    model.param().set("b_rod", b_rod + "[nm]");
    model.param().set("R_ext", Rext + "[nm]", "Radius enclosing all particles");
    model.param().set("r_sca", "1.1*R_ext + 5[nm]", "Integration sphere for sigma_sca");
    model.param().set("r_EM", "r_sca + " + g(cfg, "em_margin", "0.35") + "*lam_max/n_m", "Inner PML radius");
    model.param().set("t_pml", g(cfg, "pml_lambda", "0.4") + "*lam_max/n_m", "PML thickness");
    model.param().set("r_pml", "r_EM + t_pml");
    model.param().set("R_th", g(cfg, "R_th_factor", "100") + "*r_pml", "Radius of the thermal domain (T = T0)");
    model.param().set("h_np", hNp + "[nm]", "Mesh size in the particles");
    model.param().set("h_em", "lam_min/n_m/" + g(cfg, "el_per_lambda", "8"), "Mesh size in the medium");
    model.param().set("E0", "1[V/m]");
    model.param().set("I0", "0.5*n_m*c_const*epsilon0_const*E0^2", "Intensity of the computational field");
    model.param().set("I_phys", g(cfg, "I_phys_Wcm2", "1e5") + "[W/cm^2]", "Intensity used in the thermal problem");
    model.param().set("k_m", "2*pi*n_m/lam0");
    model.param().set("kx_p", "if(pol==3,1,0)");
    model.param().set("kz_p", "if(pol==3,0,1)");
    model.param().set("ex_p", "if(pol==1,1,0)");
    model.param().set("ey_p", "if(pol==2,1,0)");
    model.param().set("ez_p", "if(pol==3,1,0)");
    model.param().set("k_med", g(cfg, "k_med", "0.6") + "[W/(m*K)]", "Thermal conductivity of the medium");
    model.param().set("rho_med", g(cfg, "rho_med", "997") + "[kg/m^3]");
    model.param().set("cp_med", g(cfg, "cp_med", "4180") + "[J/(kg*K)]");
    model.param().set("k_np", g(cfg, "k_np", "314") + "[W/(m*K)]", "Particle thermal conductivity (default Au)");
    model.param().set("rho_np", g(cfg, "rho_np", "19300") + "[kg/m^3]");
    model.param().set("cp_np", g(cfg, "cp_np", "129") + "[J/(kg*K)]");
    model.param().set("T0", "293.15[K]");

    String matSpec = g(cfg, "material", "Au");
    String epsNP = epsExpr(model, matSpec, lamMin, lamMax, log);

    // ------------------------------------------------------------------ geometry (nm)
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("nm");
    sph(model, "s_th", "R_th");
    sph(model, "s_pml", "r_pml");
    sph(model, "s_em", "r_EM");
    sph(model, "s_sca", "r_sca");
    if (shape.equals("rod")) {
      model.component("comp1").geom("geom1").create("p0", "Ellipsoid");
      model.component("comp1").geom("geom1").feature("p0").set("semiaxes", new String[]{"a_rod", "b_rod", "b_rod"});
    } else {
      for (int k = 0; k < N; k++) {
        double[] p = P.get(k);
        model.component("comp1").geom("geom1").create("p" + k, "Sphere");
        model.component("comp1").geom("geom1").feature("p" + k).set("r", p[3] + "");
        model.component("comp1").geom("geom1").feature("p" + k).set("pos", new String[]{p[0] + "", p[1] + "", p[2] + ""});
      }
    }
    // gap refinement: auxiliary spheres (same medium) around every narrow gap, meshed with h_gap
    List<double[]> gapBalls = new ArrayList<>();
    if (cfg.containsKey("h_gap") && !shape.equals("rod")) {
      double gref = Double.parseDouble(g(cfg, "gap_ref", String.valueOf(Rmin / 2)));
      double fb = Double.parseDouble(g(cfg, "gap_ball", "1.5"));
      for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) {
        double[] a = P.get(i), b = P.get(j);
        double d = Math.sqrt(Math.pow(a[0]-b[0], 2) + Math.pow(a[1]-b[1], 2) + Math.pow(a[2]-b[2], 2));
        double gap = d - a[3] - b[3];
        if (gap <= 0 || gap > gref) continue;
        double re = 2 * a[3] * b[3] / (a[3] + b[3]);
        double rg = Math.min(0.5 * gap + fb * Math.sqrt(re * gap), 0.8 * re);
        double t = (a[3] + 0.5 * gap) / d;   // midpoint of the gap on the line of centres
        gapBalls.add(new double[]{a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1]), a[2] + t*(b[2]-a[2]), rg, gap});
      }
      for (int q = 0; q < gapBalls.size(); q++) {
        double[] gb = gapBalls.get(q);
        model.component("comp1").geom("geom1").create("gb" + q, "Sphere");
        model.component("comp1").geom("geom1").feature("gb" + q).set("r", gb[3] + "");
        model.component("comp1").geom("geom1").feature("gb" + q).set("pos", new String[]{gb[0] + "", gb[1] + "", gb[2] + ""});
        log.printf(Locale.US, "gap ball %d: centre (%.2f, %.2f, %.2f) nm, r = %.2f nm, gap = %.2f nm%n", q, gb[0], gb[1], gb[2], gb[3], gb[4]);
      }
      log.println("Gap refinement: " + gapBalls.size() + " gaps <= " + gref + " nm meshed with h_gap = " + cfg.get("h_gap") + " nm");
    }
    model.component("comp1").geom("geom1").run("fin");
    log.println("Domains: " + model.component("comp1").geom("geom1").getNDomains());

    // ------------------------------------------------------------------ selections
    String[] pnames = new String[N];
    for (int k = 0; k < N; k++) {
      double[] p = P.get(k);
      double rb = shape.equals("rod") ? a_rod * 1.02 : p[3] * (gapBalls.isEmpty() ? 1.02 : 1.002);
      ball(model, "sb" + k, 3, p[0], p[1], p[2], rb + "", "Particle " + k);
      pnames[k] = "sb" + k;
    }
    model.component("comp1").selection().create("sel_np", "Union");
    model.component("comp1").selection("sel_np").set("entitydim", 3);
    model.component("comp1").selection("sel_np").set("input", pnames);
    ball(model, "sel_inner", 3, 0, 0, 0, "r_EM*1.001", "Inside r_EM");
    ball(model, "sel_emw", 3, 0, 0, 0, "r_pml*1.001", "EM domain (with PML)");
    model.component("comp1").selection().create("sel_pml", "Difference");
    model.component("comp1").selection("sel_pml").set("entitydim", 3);
    model.component("comp1").selection("sel_pml").set("add", new String[]{"sel_emw"});
    model.component("comp1").selection("sel_pml").set("subtract", new String[]{"sel_inner"});
    model.component("comp1").selection().create("sel_med_in", "Difference");
    model.component("comp1").selection("sel_med_in").set("entitydim", 3);
    model.component("comp1").selection("sel_med_in").set("add", new String[]{"sel_inner"});
    model.component("comp1").selection("sel_med_in").set("subtract", new String[]{"sel_np"});
    ball(model, "sel_sca_out", 2, 0, 0, 0, "r_sca*1.001", "");
    ball(model, "sel_sca_in", 2, 0, 0, 0, "r_sca*0.999", "");
    model.component("comp1").selection().create("sel_sca", "Difference");
    model.component("comp1").selection("sel_sca").label("Integration sphere for sigma_sca");
    model.component("comp1").selection("sel_sca").set("entitydim", 2);
    model.component("comp1").selection("sel_sca").set("add", new String[]{"sel_sca_out"});
    model.component("comp1").selection("sel_sca").set("subtract", new String[]{"sel_sca_in"});
    model.component("comp1").selection().create("sel_all", "Explicit");
    model.component("comp1").selection("sel_all").all();
    model.component("comp1").selection().create("sel_ext", "Adjacent");
    model.component("comp1").selection("sel_ext").set("input", new String[]{"sel_all"});
    model.component("comp1").selection("sel_ext").set("outputdim", 2);
    model.component("comp1").selection("sel_ext").set("interior", false);
    model.component("comp1").selection("sel_ext").set("exterior", true);
    try {
      String[] ss = {"sel_np", "sel_inner", "sel_pml", "sel_med_in", "sel_sca", "sel_ext"};
      int[] dd = {3, 3, 3, 3, 2, 2};
      for (int q = 0; q < ss.length; q++) log.println("Selection " + ss[q] + ": " + Arrays.toString(model.component("comp1").selection(ss[q]).entities(dd[q])));
    } catch (Exception e) { log.println("WARN diag: " + e.getMessage()); }

    // ------------------------------------------------------------------ variables and operators
    model.component("comp1").variable().create("var_g");
    model.component("comp1").variable("var_g").set("Ebx", "E0*ex_p*exp(-i*k_m*(kx_p*x + kz_p*z))");
    model.component("comp1").variable("var_g").set("Eby", "E0*ey_p*exp(-i*k_m*(kx_p*x + kz_p*z))");
    model.component("comp1").variable("var_g").set("Ebz", "E0*ez_p*exp(-i*k_m*(kx_p*x + kz_p*z))");
    model.component("comp1").variable("var_g").set("E2n", "abs(emw.Ex)^2 + abs(emw.Ey)^2 + abs(emw.Ez)^2");
    model.component("comp1").variable("var_g").set("sig_abs", "intop_np(emw.Qh)/I0");
    model.component("comp1").variable("var_g").set("sig_sca",
        "intop_sca(0.5*real((emw.relEy*conj(emw.relHz)-emw.relEz*conj(emw.relHy))*x"
      + " + (emw.relEz*conj(emw.relHx)-emw.relEx*conj(emw.relHz))*y"
      + " + (emw.relEx*conj(emw.relHy)-emw.relEy*conj(emw.relHx))*z)/sqrt(x^2+y^2+z^2))/I0");
    model.component("comp1").variable("var_g").set("P_abs_phys", "intop_np(emw.Qh)*I_phys/I0", "Absorbed power at I_phys");
    // explicit Joule/polarization loss density, valid in the thermal study where emw.freq is undefined
    model.component("comp1").variable().create("var_np");
    model.component("comp1").variable("var_np").selection().geom("geom1", 3);
    model.component("comp1").variable("var_np").selection().named("sel_np");
    model.component("comp1").variable("var_np").set("Q_np",
        "0.5*(2*pi*c_const/lam0)*epsilon0_const*" + imagEps(epsNP) + "*(abs(emw.Ex)^2 + abs(emw.Ey)^2 + abs(emw.Ez)^2)",
        "Heat source density in the particles (= emw.Qh)");
    intop(model, "intop_np", "sel_np");
    intop(model, "intop_sca", "sel_sca");
    for (int k = 0; k < N; k++) {
      intop(model, "ip" + k, "sb" + k);
      model.component("comp1").cpl().create("av" + k, "Average");
      model.component("comp1").cpl("av" + k).set("opname", "av" + k);
      model.component("comp1").cpl("av" + k).selection().named("sb" + k);
    }
    model.component("comp1").cpl().create("max_np", "Maximum");
    model.component("comp1").cpl("max_np").set("opname", "max_np");
    model.component("comp1").cpl("max_np").selection().named("sel_np");

    // ------------------------------------------------------------------ materials
    String one9 = "1";
    model.component("comp1").material().create("mat_med", "Common");
    model.component("comp1").material("mat_med").label("Medium");
    model.component("comp1").material("mat_med").selection().all();
    iso(model, "mat_med", "relpermittivity", "n_m^2");
    iso(model, "mat_med", "relpermeability", one9);
    iso(model, "mat_med", "electricconductivity", "0");
    iso(model, "mat_med", "thermalconductivity", "k_med");
    model.component("comp1").material("mat_med").propertyGroup("def").set("density", "rho_med");
    model.component("comp1").material("mat_med").propertyGroup("def").set("heatcapacity", "cp_med");
    model.component("comp1").material().create("mat_au", "Common");
    model.component("comp1").material("mat_au").label("Particles: " + matSpec);
    model.component("comp1").material("mat_au").selection().named("sel_np");
    iso(model, "mat_au", "relpermittivity", epsNP);
    iso(model, "mat_au", "relpermeability", one9);
    iso(model, "mat_au", "electricconductivity", "0");
    iso(model, "mat_au", "thermalconductivity", "k_np");
    model.component("comp1").material("mat_au").propertyGroup("def").set("density", "rho_np");
    model.component("comp1").material("mat_au").propertyGroup("def").set("heatcapacity", "cp_np");

    // ------------------------------------------------------------------ PML
    model.component("comp1").coordSystem().create("pml1", "PML");
    model.component("comp1").coordSystem("pml1").selection().named("sel_pml");
    try { model.component("comp1").coordSystem("pml1").set("ScalingType", "Spherical"); } catch (Exception e) { log.println("WARN PML type: " + e.getMessage()); }

    // ------------------------------------------------------------------ EM physics
    model.component("comp1").physics().create("ewfd", "ElectromagneticWaves", "geom1");
    model.component("comp1").physics("ewfd").selection().named("sel_emw");
    try { model.component("comp1").physics("ewfd").prop("BackgroundField").set("SolveFor", "scatteredField"); }
    catch (Exception e) { log.println("WARN SolveFor: " + e.getMessage()); }
    try { model.component("comp1").physics("ewfd").prop("BackgroundField").set("Eb", new String[]{"Ebx", "Eby", "Ebz"}); }
    catch (Exception e) { log.println("WARN Eb: " + e.getMessage()); }

    // ------------------------------------------------------------------ heat transfer
    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T0");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 3);
    model.component("comp1").physics("ht").feature("hs1").label("Absorption (scaled emw.Qh)");
    model.component("comp1").physics("ht").feature("hs1").selection().named("sel_np");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q_np*I_phys/I0");
    model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 2);
    model.component("comp1").physics("ht").feature("temp1").selection().named("sel_ext");
    model.component("comp1").physics("ht").feature("temp1").set("T0", "T0");

    // ------------------------------------------------------------------ mesh
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "R_th/5");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "h_np/4");
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", g(cfg, "hgrad", "1.4"));
    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    size(model, "sz_np", 3, "sel_np", "h_np");
    if (cfg.containsKey("h_surf")) {   // fine mesh only on the particle surfaces (gaps), coarser inside
      model.param().set("h_surf", cfg.get("h_surf") + "[nm]", "Mesh size on the particle surfaces");
      model.component("comp1").selection().create("sel_np_b", "Adjacent");
      model.component("comp1").selection("sel_np_b").set("input", new String[]{"sel_np"});
      model.component("comp1").selection("sel_np_b").set("outputdim", 2);
      size(model, "sz_surf", 2, "sel_np_b", "h_surf");
    }
    size(model, "sz_in", 3, "sel_med_in", "h_em");
    size(model, "sz_pml", 3, "sel_pml", "h_em");
    if (!gapBalls.isEmpty()) {   // last Size node wins on the gap domains (medium and particle caps)
      model.param().set("h_gap", cfg.get("h_gap") + "[nm]", "Mesh size in the gap regions");
      String[] gn = new String[gapBalls.size()];
      for (int q = 0; q < gapBalls.size(); q++) {
        double[] gb = gapBalls.get(q);
        ball(model, "sg" + q, 3, gb[0], gb[1], gb[2], (gb[3] * 1.002) + "", "");
        gn[q] = "sg" + q;
      }
      model.component("comp1").selection().create("sel_gap", "Union");
      model.component("comp1").selection("sel_gap").set("entitydim", 3);
      model.component("comp1").selection("sel_gap").set("input", gn);
      size(model, "sz_gap", 3, "sel_gap", "h_gap");
      try { log.println("Selection sel_gap: " + Arrays.toString(model.component("comp1").selection("sel_gap").entities(3))); } catch (Exception e) { }
    }
    model.component("comp1").mesh("mesh1").run();
    try { log.println("Mesh elements: " + model.component("comp1").mesh("mesh1").getNumElem()); } catch (Exception e) { }

    // ------------------------------------------------------------------ studies
    model.study().create("std1");
    model.study("std1").create("freq", "Frequency");
    model.study("std1").feature("freq").set("plist", "c_const/lam0");
    model.study("std1").feature("freq").setEntry("activate", "ht", false);
    model.sol().create("sol1"); model.sol("sol1").study("std1"); model.sol("sol1").attach("std1"); model.sol("sol1").createAutoSequence("std1");
    // thermal study = two steps in ONE study (frequency domain -> stationary heat): the second step
    // inherits the full EM solution. A separate stationary study pointing to std1 only received the
    // background field (heat source ~1/|E_in/E0|^2 too small), see AGENTS.md.
    model.study().create("std2");
    model.study("std2").create("freq", "Frequency");
    model.study("std2").feature("freq").set("plist", "c_const/lam0");
    model.study("std2").feature("freq").setEntry("activate", "ht", false);
    model.study("std2").create("stat", "Stationary");
    model.study("std2").feature("stat").setEntry("activate", "ewfd", false);

    model.result().dataset().create("dset_em", "Solution");
    model.result().dataset("dset_em").set("solution", "sol1");
    model.result().numerical().create("gev", "EvalGlobal");
    model.result().numerical("gev").set("data", "dset_em");

    // ---- 1) spectra
    PrintWriter sp = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_spectra.txt"));
    sp.println("% lambda_nm  pol  sig_abs_nm2  sig_sca_nm2  sig_ext_nm2");
    for (double pd : pols) {
      int p = (int) pd;
      for (double L : lams) {
        model.param().set("pol", String.valueOf(p));
        model.param().set("lam0", L + "[nm]");
        long t0 = System.currentTimeMillis();
        model.sol("sol1").runAll();
        double[] v = evalSafe(model, new String[]{"sig_abs", "sig_sca"}, new String[]{"nm^2", "nm^2"}, log);
        sp.printf(Locale.US, "%8.2f %3d %14.6e %14.6e %14.6e%n", L, p, v[0], v[1], v[0] + v[1]);
        sp.flush();
        log.printf(Locale.US, "lam %.1f pol %d: sig_abs %.4e sig_sca %.4e nm^2 (%.0f s)%n", L, p, v[0], v[1], (System.currentTimeMillis() - t0) / 1e3);
      }
    }
    sp.close();

    // ---- 2) field and heat at lam_heat, pol_heat
    model.param().set("pol", String.valueOf(polHeat));
    model.param().set("lam0", lamHeat + "[nm]");
    model.sol("sol1").runAll();
    double rEM = model.param().evaluate("r_EM") * 1e9;
    model.result().dataset().create("cl_x", "CutLine3D");
    model.result().dataset("cl_x").set("data", "dset_em");
    model.result().dataset("cl_x").set("genpoints", new double[][]{{-rEM, 0, 0}, {rEM, 0, 0}});
    model.result().export().create("ex_E", "Data");
    model.result().export("ex_E").set("data", "cl_x");
    model.result().export("ex_E").set("expr", new String[]{"E2n/E0^2"});
    model.result().export("ex_E").set("unit", new String[]{"1"});
    model.result().export("ex_E").set("filename", BASE + "results/" + tag + "_E2_xline.txt");
    model.result().export("ex_E").run();

    sol(model, "sol2", "std2");
    model.result().dataset().create("dset_ht", "Solution");
    model.result().dataset("dset_ht").set("solution", "sol2");
    model.result().numerical("gev").set("data", "dset_em");
    double[] em = evalSafe(model, new String[]{"P_abs_phys", "sig_abs", "intop_np(Q_np)/intop_np(emw.Qh)", "intop_np(1)"}, new String[]{"W", "nm^2", "1", "nm^3"}, log);
    log.printf(Locale.US, "check: int(Q_np)/int(emw.Qh) = %.6f (must be 1)%n", em[2]);
    double vAn = 0; for (double[] p : P) vAn += 4.0 / 3 * Math.PI * p[3] * p[3] * p[3];
    if (shape.equals("rod")) vAn = 4.0 / 3 * Math.PI * a_rod * b_rod * b_rod;
    log.printf(Locale.US, "check: particle volume FEM/analytic = %.4f (selection of the metal, must be ~1)%n", em[3] / vAn);
    String[] dx = {"intop_np(abs(Ex)^2)/intop_np(abs(emw.Ex)^2)", "intop_np(abs(emw.relEx)^2)/intop_np(abs(emw.Ex)^2)"};
    for (String q : dx) log.printf(Locale.US, "diag %s = %.6f%n", q, evalSafe(model, new String[]{q}, new String[]{"1"}, log)[0]);
    model.result().numerical("gev").set("data", "dset_ht");
    double[] tt = evalSafe(model, new String[]{"max_np(T)-T0", "intop_np(Q_np)*I_phys/I0"}, new String[]{"K", "W"}, log);
    log.printf(Locale.US, "check: heat source seen by the thermal study = %.6e W (EM study: %.6e W)%n", tt[1], em[0]);
    if (!(tt[1] > 0.5 * em[0])) log.println("WARN the thermal study does not see the EM heat source");
    double[] th = {em[0], tt[0], em[1]};
    double Rth = model.param().evaluate("R_th");
    PrintWriter tf = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_thermal.txt"));
    tf.println("% lambda_nm pol I_phys_Wcm2 P_abs_W dT_max_K sig_abs_nm2 R_th_m k_med");
    tf.printf(Locale.US, "%8.2f %3d %12.4e %14.6e %14.6e %14.6e %12.4e %8.4f%n", lamHeat, polHeat,
        model.param().evaluate("I_phys") / 1e4, th[0], th[1], th[2], Rth, model.param().evaluate("k_med"));
    tf.close();
    // per particle: x y z R sig_abs_i dT_mean_i
    PrintWriter pf = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_particles.txt"));
    pf.println("% i x_nm y_nm z_nm R_nm sig_abs_i_nm2 dT_mean_K   (lambda = " + lamHeat + " nm, pol = " + polHeat + ")");
    for (int k = 0; k < N; k++) {
      double[] p = P.get(k);
      model.result().numerical("gev").set("data", "dset_em");
      double[] v1 = evalSafe(model, new String[]{"ip" + k + "(emw.Qh)/I0"}, new String[]{"nm^2"}, log);
      model.result().numerical("gev").set("data", "dset_ht");
      double[] v2 = evalSafe(model, new String[]{"av" + k + "(T)-T0"}, new String[]{"K"}, log);
      double[] v = {v1[0], v2[0]};
      pf.printf(Locale.US, "%4d %10.3f %10.3f %10.3f %8.3f %14.6e %14.6e%n", k, p[0], p[1], p[2], p[3], v[0], v[1]);
    }
    pf.close();
    model.result().dataset().create("cl_T", "CutLine3D");
    model.result().dataset("cl_T").set("data", "dset_ht");
    model.result().dataset("cl_T").set("genpoints", new double[][]{{0, 0, 0}, {50 * Rext, 0, 0}});
    model.result().export().create("ex_T", "Data");
    model.result().export("ex_T").set("data", "cl_T");
    model.result().export("ex_T").set("expr", new String[]{"T-T0"});
    model.result().export("ex_T").set("unit", new String[]{"K"});
    model.result().export("ex_T").set("filename", BASE + "results/" + tag + "_T_xline.txt");
    model.result().export("ex_T").run();

    // ---- plots for the GUI
    model.result().create("pg_E", "PlotGroup3D");
    model.result("pg_E").label("|E|^2/E0^2");
    model.result("pg_E").set("data", "dset_em");
    model.result("pg_E").create("slc1", "Slice");
    model.result("pg_E").feature("slc1").set("expr", "E2n/E0^2");
    model.result("pg_E").feature("slc1").set("quickplane", "xy");
    model.result("pg_E").feature("slc1").set("quickznumber", 1);
    model.result().create("pg_T", "PlotGroup3D");
    model.result("pg_T").label("DeltaT");
    model.result("pg_T").set("data", "dset_ht");
    model.result("pg_T").create("slc1", "Slice");
    model.result("pg_T").feature("slc1").set("expr", "T-T0");
    model.result("pg_T").feature("slc1").set("quickplane", "xy");
    model.result("pg_T").feature("slc1").set("quickznumber", 1);
    log.println("END");
    log.close();
    return model;
  }

  // ---------------------------------------------------------------- helper
  static void sph(Model m, String tag, String r) {
    m.component("comp1").geom("geom1").create(tag, "Sphere");
    m.component("comp1").geom("geom1").feature(tag).set("r", r);
  }
  static void ball(Model m, String tag, int dim, double x, double y, double z, String r, String label) {
    m.component("comp1").selection().create(tag, "Ball");
    if (!label.isEmpty()) m.component("comp1").selection(tag).label(label);
    m.component("comp1").selection(tag).set("entitydim", dim);
    m.component("comp1").selection(tag).set("posx", x + "");
    m.component("comp1").selection(tag).set("posy", y + "");
    m.component("comp1").selection(tag).set("posz", z + "");
    m.component("comp1").selection(tag).set("r", r);
    m.component("comp1").selection(tag).set("condition", "inside");
  }
  static void intop(Model m, String tag, String sel) {
    m.component("comp1").cpl().create(tag, "Integration");
    m.component("comp1").cpl(tag).set("opname", tag);
    m.component("comp1").cpl(tag).selection().named(sel);
  }
  static void iso(Model m, String mat, String prop, String v) {
    m.component("comp1").material(mat).propertyGroup("def").set(prop, new String[]{v, "0", "0", "0", v, "0", "0", "0", v});
  }
  static void size(Model m, String tag, int dim, String sel, String hmax) {
    m.component("comp1").mesh("mesh1").feature("ftet1").create(tag, "Size");
    m.component("comp1").mesh("mesh1").feature("ftet1").feature(tag).selection().geom("geom1", dim);
    m.component("comp1").mesh("mesh1").feature("ftet1").feature(tag).selection().named(sel);
    m.component("comp1").mesh("mesh1").feature("ftet1").feature(tag).set("custom", "on");
    m.component("comp1").mesh("mesh1").feature("ftet1").feature(tag).set("hmaxactive", true);
    m.component("comp1").mesh("mesh1").feature("ftet1").feature(tag).set("hmax", hmax);
  }
  static void sol(Model m, String s, String std) {
    m.sol().create(s); m.sol(s).study(std); m.sol(s).attach(std); m.sol(s).createAutoSequence(std); m.sol(s).runAll();
  }
  static double[] evalSafe(Model m, String[] ex, String[] un, PrintWriter log) {
    double[] out = new double[ex.length];
    for (int k = 0; k < ex.length; k++) {
      try {
        m.result().numerical("gev").set("expr", new String[]{ex[k]});
        m.result().numerical("gev").set("unit", new String[]{un[k]});
        out[k] = m.result().numerical("gev").getReal()[0][0];
      } catch (Exception e) { out[k] = Double.NaN; log.println("WARN eval " + ex[k] + ": " + e.getMessage()); }
    }
    return out;
  }

  // ---------------------------------------------------------------- materials (shared by the templates)
  // spec: "Au" (analytic, Etchegoin 2006) | "nk:<file relative to the repository>" (lambda_nm n k) | a real index "1.45"
  static int nkCount = 0;
  static boolean auDone = false;
  static double[] lastNK = {0, 0};
  static String epsExpr(Model model, String spec, double lamMin, double lamMax, PrintWriter log) throws IOException {
    spec = spec.trim();
    if (spec.equalsIgnoreCase("Au")) {
      if (!auDone) {
        String et = "(1.53 - 1/(145^2*(1/lam^2 + i/(17000*lam)))"
          + " + 0.94/468*(exp(-i*pi/4)/(1/468 - 1/lam - i/2300) + exp(i*pi/4)/(1/468 + 1/lam + i/2300))"
          + " + 1.36/331*(exp(-i*pi/4)/(1/331 - 1/lam - i/940) + exp(i*pi/4)/(1/331 + 1/lam + i/940)))";
        model.func().create("au_re", "Analytic");
        model.func("au_re").set("funcname", "epsAu_re");
        model.func("au_re").set("args", new String[]{"lam"});
        model.func("au_re").set("expr", "real" + et);
        model.func().create("au_im", "Analytic");
        model.func("au_im").set("funcname", "epsAu_im");
        model.func("au_im").set("args", new String[]{"lam"});
        model.func("au_im").set("expr", "imag" + et);
        auDone = true;
      }
      lastNK = new double[]{0.5, 2.0};
      return "epsAu_re(lam_nm) - i*epsAu_im(lam_nm)";
    }
    if (spec.startsWith("nk:")) {
      String path = BASE + spec.substring(3).trim();
      List<String[]> rows = new ArrayList<>();
      BufferedReader br = new BufferedReader(new FileReader(path));
      String s;
      while ((s = br.readLine()) != null) {
        s = s.trim(); if (s.isEmpty() || s.startsWith("#") || s.startsWith("%")) continue;
        String[] t = s.split("[\\s,;]+"); if (t.length < 3) continue;
        try { Double.parseDouble(t[0]); } catch (NumberFormatException e) { continue; }
        rows.add(new String[]{t[0], t[1], t[2]});
      }
      br.close();
      double l0 = Double.parseDouble(rows.get(0)[0]), l1 = Double.parseDouble(rows.get(rows.size() - 1)[0]);
      if (lamMin < Math.min(l0, l1) || lamMax > Math.max(l0, l1))
        log.println("WARN " + spec + ": table range " + l0 + "-" + l1 + " nm does not cover " + lamMin + "-" + lamMax + " nm (constant extrapolation)");
      String tag = "nk" + (nkCount++);
      String[][] tn = new String[rows.size()][2], tk = new String[rows.size()][2];
      double best = 1e9;
      for (int q = 0; q < rows.size(); q++) {
        tn[q][0] = rows.get(q)[0]; tn[q][1] = rows.get(q)[1];
        tk[q][0] = rows.get(q)[0]; tk[q][1] = rows.get(q)[2];
        double d = Math.abs(Double.parseDouble(rows.get(q)[0]) - lamMin);
        if (d < best) { best = d; lastNK = new double[]{Double.parseDouble(rows.get(q)[1]), Double.parseDouble(rows.get(q)[2])}; }
      }
      for (String[] fk : new String[][]{{tag + "_n", "n"}, {tag + "_k", "k"}}) {
        model.func().create(fk[0], "Interpolation");
        model.func(fk[0]).label(spec + " (" + fk[1] + ")");
        model.func(fk[0]).set("funcname", fk[0]);
        model.func(fk[0]).set("table", fk[1].equals("n") ? tn : tk);
        try { model.func(fk[0]).set("interp", "linear"); } catch (Exception e) { log.println("WARN interp: " + e.getMessage()); }   // same as Python
        try { model.func(fk[0]).set("extrap", "const"); } catch (Exception e) { log.println("WARN extrap: " + e.getMessage()); }
      }
      log.println("material " + spec + ": " + rows.size() + " points, " + l0 + "-" + l1 + " nm -> functions " + tag + "_n, " + tag + "_k");
      return "(" + tag + "_n(lam_nm) - i*" + tag + "_k(lam_nm))^2";
    }
    double n = Double.parseDouble(spec);
    lastNK = new double[]{n, 0};
    return "(" + spec + ")^2";
  }
  static String imagEps(String epsExpr) {   // Im part with the exp(-iwt) sign (positive for losses)
    return "(-imag(" + epsExpr + "))";
  }

  public static void main(String[] args) { run(); }
}
