/*
 * PT_Multilayer.java  -  ptpipe COMSOL template, level L2 (planar stacks / unit cell, nm)
 *
 * 2D unit cell of width w_cell, Floquet periodic conditions on the lateral edges and
 * periodic ports on top (excited) and bottom, as in the numerical ellipsometry approach of
 * Lio et al., RSC Adv. 9, 21429 (2019). For a laterally uniform stack the result must coincide
 * with the transfer-matrix method (ptpipe.optics.tmm); a patterned layer (metasurface) can be
 * added to the same cell.
 *
 * Output per wavelength: R = |S11|^2, T = |S21|^2, absorption in every layer from the
 * dissipated power int(Qh) (independent of R and T: energy conservation is a check), and the
 * heat-source profile a(z) = q(z)/I_inc along the cell axis at lam_profile.
 *
 * Configuration: comsol/multilayer.cfg
 *   layers = Au:15, 1.45:100, nk:materials/Ti.txt:5   (top -> bottom; material:thickness_nm,
 *            material = Au | nk:<file> | real index)
 *   n_sup, n_sub, lams (a:step:b or list), lam_profile, theta_deg, pol (s | p), w_cell
 * Output: results/<tag>_RTA.txt, results/<tag>_qz.txt, results/<tag>_log.txt
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class PT_Multilayer {

  static String BASE = "/path/to/ptpipe/";   // set automatically by run_validation.sh

  public static Model run() {
    try { return runImpl(); } catch (Exception e) { throw new RuntimeException(e); }
  }

  static Model runImpl() throws Exception {
    try { String r = System.getProperty("ptpipe.root"); if (r != null) BASE = r.endsWith("/") ? r : r + "/"; } catch (Exception e) { }
    Map<String, String> cfg = readCfg(BASE + "comsol/multilayer.cfg");
    String tag = g(cfg, "tag", "stack");
    PrintWriter log = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_log.txt"), true);
    log.println("cfg: " + cfg);
    double[] lams = list(g(cfg, "lams", "450:10:900"));
    double lamProf = Double.parseDouble(g(cfg, "lam_profile", "540"));
    double lamMin = Math.min(Arrays.stream(lams).min().getAsDouble(), lamProf);
    double lamMax = Math.max(Arrays.stream(lams).max().getAsDouble(), lamProf);
    String pol = g(cfg, "pol", "s");

    // layers top -> bottom
    String[] toks = g(cfg, "layers", "Au:15, 1.45:100, Au:100").split(",");
    int NL = toks.length;
    String[] mat = new String[NL]; double[] d = new double[NL];
    for (int k = 0; k < NL; k++) {          // "<material spec>:<thickness nm>", spec may itself contain ':' (nk:file)
      String t = toks[k].trim(); int c = t.lastIndexOf(':');
      mat[k] = t.substring(0, c).trim(); d[k] = Double.parseDouble(t.substring(c + 1).trim());
    }
    double D = 0; for (double x : d) D += x;

    Model model = ModelUtil.create("Model");
    model.modelPath(BASE + "comsol");
    model.label("PT_Multilayer_" + tag + ".mph");
    model.param().set("lam0", lamProf + "[nm]");
    model.param().set("lam_nm", "lam0/1[nm]");
    model.param().set("lam_min", lamMin + "[nm]");
    model.param().set("lam_max", lamMax + "[nm]");
    model.param().set("n_sup", g(cfg, "n_sup", "1.0"));
    model.param().set("n_sub", g(cfg, "n_sub", "1.52"));
    model.param().set("theta", g(cfg, "theta_deg", "0") + "[deg]", "Angle of incidence (in the superstrate)");
    model.param().set("theta_t", g(cfg, "port2_sign", "1") + "*asin(n_sup*sin(theta)/n_sub)",
        "Angle of the transmitted wave in the substrate (Snell); the listener port must use it");
    model.param().set("w_cell", g(cfg, "w_cell", "50") + "[nm]", "Cell width (sub-wavelength: only the 0th order)");
    model.param().set("h_sup", "0.6*lam_max/n_sup", "Superstrate height above the stack");
    model.param().set("h_sub", "0.6*lam_max/n_sub", "Substrate height below the stack");
    model.param().set("D_stack", D + "[nm]");
    model.param().set("P_in", "1[W]", "Port input power (per 1 m out of plane)");

    // ------------------------------------------------------------------ geometry
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 2);
    model.component("comp1").geom("geom1").lengthUnit("nm");
    rect(model, "r_sub", "-h_sub", "h_sub");
    double y = D;
    double[] y0 = new double[NL], y1 = new double[NL];
    for (int k = 0; k < NL; k++) { y1[k] = y; y0[k] = y - d[k]; y = y0[k]; rect(model, "r_L" + k, y0[k] + "", d[k] + ""); }
    rect(model, "r_sup", "D_stack", "h_sup");
    model.component("comp1").geom("geom1").run("fin");

    // ------------------------------------------------------------------ selections (boxes, robust to numbering)
    double eps = 1e-3;
    double hsup = model.param().evaluate("h_sup") * 1e9, hsub = model.param().evaluate("h_sub") * 1e9;
    double w = model.param().evaluate("w_cell") * 1e9;
    String[] Lsel = new String[NL];
    for (int k = 0; k < NL; k++) { box(model, "sL" + k, 2, -eps, w + eps, y0[k] - eps, y1[k] + eps); Lsel[k] = "sL" + k; }
    box(model, "s_sup", 2, -eps, w + eps, D - eps, D + hsup + eps);
    box(model, "s_sub", 2, -eps, w + eps, -hsub - eps, eps);
    box(model, "b_top", 1, -eps, w + eps, D + hsup - eps, D + hsup + eps);
    box(model, "b_bot", 1, -eps, w + eps, -hsub - eps, -hsub + eps);
    box(model, "b_left", 1, -eps, eps, -hsub - eps, D + hsup + eps);
    box(model, "b_right", 1, w - eps, w + eps, -hsub - eps, D + hsup + eps);
    model.component("comp1").selection().create("b_lr", "Union");
    model.component("comp1").selection("b_lr").set("entitydim", 1);
    model.component("comp1").selection("b_lr").set("input", new String[]{"b_left", "b_right"});
    box(model, "b_all", 1, -eps, w + eps, -hsub - eps, D + hsup + eps);
    model.component("comp1").selection().create("b_hor", "Difference");
    model.component("comp1").selection("b_hor").set("entitydim", 1);
    model.component("comp1").selection("b_hor").set("add", new String[]{"b_all"});
    model.component("comp1").selection("b_hor").set("subtract", new String[]{"b_lr"});
    model.component("comp1").selection().create("s_stack", "Union");
    model.component("comp1").selection("s_stack").set("entitydim", 2);
    model.component("comp1").selection("s_stack").set("input", Lsel);

    // ------------------------------------------------------------------ materials
    material(model, "m_sup", "s_sup", "n_sup^2", null);
    material(model, "m_sub", "s_sub", "n_sub^2", null);
    double[] nkLay = new double[NL], kLay = new double[NL];
    for (int k = 0; k < NL; k++) {
      String e = epsExpr(model, mat[k], lamMin, lamMax, log);
      nkLay[k] = lastNK[0]; kLay[k] = lastNK[1];
      material(model, "m_L" + k, "sL" + k, e, null);
      model.component("comp1").material("m_L" + k).label("Layer " + k + ": " + mat[k] + " " + d[k] + " nm");
    }

    // ------------------------------------------------------------------ variables
    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Rr", "abs(emw.S11)^2", "Reflectance");
    model.component("comp1").variable("var1").set("Tt", "abs(emw.S21)^2", "Transmittance");
    model.component("comp1").variable("var1").set("I_inc", "P_in/(w_cell*1[m])", "Incident intensity (normal incidence)");
    for (int k = 0; k < NL; k++) {
      intop(model, "iL" + k, "sL" + k);
      model.component("comp1").variable("var1").set("A_L" + k, "iL" + k + "(emw.Qh)*1[m]/P_in", "Absorption in layer " + k);
    }
    intop(model, "i_stack", "s_stack");
    model.component("comp1").variable("var1").set("A_Q", "i_stack(emw.Qh)*1[m]/P_in", "Total absorption from Qh");

    // ------------------------------------------------------------------ physics
    model.component("comp1").physics().create("ewfd", "ElectromagneticWaves", "geom1");
    model.component("comp1").physics("ewfd").create("port1", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port1").selection().named("b_top");
    model.component("comp1").physics("ewfd").create("port2", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port2").selection().named("b_bot");
    for (String pt : new String[]{"port1", "port2"}) {
      try { model.component("comp1").physics("ewfd").feature(pt).set("PortType", "Periodic"); } catch (Exception e) { log.println("WARN PortType " + pt + ": " + e.getMessage()); }
    }
    try { model.component("comp1").physics("ewfd").feature("port1").set("PortExcitation", "on"); } catch (Exception e) { log.println("WARN PortExcitation: " + e.getMessage()); }
    try { model.component("comp1").physics("ewfd").feature("port1").set("Pin", "P_in"); } catch (Exception e) { log.println("WARN Pin: " + e.getMessage()); }
    // angle of incidence: in 2D COMSOL reads alpha_inc (alpha1_inc/alpha2_inc are the 3D elevation/azimuth; the 6.4
    // property dump of mim_s45 showed alpha1_inc = theta but alpha_inc = 0, i.e. normal incidence)
    try { model.component("comp1").physics("ewfd").feature("port1").set("alpha_inc", "theta"); } catch (Exception e) { log.println("WARN alpha_inc: " + e.getMessage()); }
    // the listener port in the substrate carries the refracted wave: same tangential k, angle theta_t
    // (with alpha_inc = 0 there, R and T at 45 deg matched no TMM angle: max |dR| = 0.19)
    try { model.component("comp1").physics("ewfd").feature("port2").set("alpha_inc", "theta_t"); } catch (Exception e) { log.println("WARN port2 alpha_inc: " + e.getMessage()); }
    // mode field of the periodic ports: s = E out of plane; p = E in plane, orthogonal to k
    // mode field of the periodic ports. In COMSOL 6.4 the user-defined polarization is read from Eampl
    // (default {t1x, t1y, t1z}: the port tangent, i.e. in-plane), not from E0: with E0 only, s and p gave
    // identical spectra at 45 deg, and normal incidence passed by coincidence. For p the in-plane field
    // must be orthogonal to k: (cos a, -sin a, 0) for alpha_inc = a with this geometry (verified against
    // the TMM at 45 deg: |dR| = 2e-3; the opposite sign gives |dR| = 0.95). The angle sign of port2 does
    // not change the result.
    String[] pts = {"port1", "port2"}, ang = {"theta", "theta_t"};
    String psgn = g(cfg, "p_sign", "-1");
    for (int q = 0; q < 2; q++) {
      String[] E0 = pol.equals("s") ? new String[]{"0", "0", "1"}
                                    : new String[]{"cos(" + ang[q] + ")", psgn + "*sin(" + ang[q] + ")", "0"};
      try {
        model.component("comp1").physics("ewfd").feature(pts[q]).set("Polarization", "UserDefined");
        model.component("comp1").physics("ewfd").feature(pts[q]).set("Eampl", E0);
        model.component("comp1").physics("ewfd").feature(pts[q]).set("E0", E0);
      } catch (Exception e) { log.println("WARN polarization " + pts[q] + ": " + e.getMessage()); }
      for (String pr : new String[]{"Polarization", "LinearPol", "Eampl", "E0", "alpha_inc"}) {
        try { log.println(pts[q] + " " + pr + " = " + String.join(", ", model.component("comp1").physics("ewfd").feature(pts[q]).getStringArray(pr))); }
        catch (Exception e) { try { log.println(pts[q] + " " + pr + " = " + model.component("comp1").physics("ewfd").feature(pts[q]).getString(pr)); } catch (Exception e2) { } }
      }
    }
    // diagnostic: list the settings of the excited port (names differ between COMSOL versions)
    try {
      String[] pr = model.component("comp1").physics("ewfd").feature("port1").properties();
      for (String q : pr) {
        String v;
        try { v = model.component("comp1").physics("ewfd").feature("port1").getString(q); } catch (Exception e) { v = "(array)"; }
        log.println("port1 property " + q + " = " + v);
      }
    } catch (Exception e) { log.println("WARN port property dump: " + e.getMessage()); }
    try { log.println("ewfd components = " + model.component("comp1").physics("ewfd").prop("components").getString("components")); } catch (Exception e) { log.println("ewfd components: n/a (" + e.getMessage() + ")"); }
    model.component("comp1").physics("ewfd").create("pc1", "PeriodicCondition", 1);
    model.component("comp1").physics("ewfd").feature("pc1").selection().named("b_lr");
    try { model.component("comp1").physics("ewfd").feature("pc1").set("PeriodicType", "Floquet"); } catch (Exception e) { log.println("WARN PeriodicType: " + e.getMessage()); }
    try { model.component("comp1").physics("ewfd").feature("pc1").set("Floquet_source", "FromPeriodicPort"); } catch (Exception e) { log.println("WARN Floquet_source: " + e.getMessage()); }

    // ------------------------------------------------------------------ mapped mesh (identical on the periodic edges)
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").create("map1", "Map");
    distBox(model, "dh", "b_hor", 4);
    for (int k = 0; k < NL; k++) {
      // lossy (k > 0.5): resolve the skin depth; dielectric: lambda/(16 n)
      double h = kLay[k] > 0.5 ? Math.min(1.5, lamMin / (4 * Math.PI * kLay[k]) / 8) : lamMin / Math.max(nkLay[k], 1.0) / 16;
      int n = Math.max(4, (int) Math.ceil(d[k] / h));
      vbox(model, "eL" + k, w, y0[k], y1[k]);   // only the two vertical edges of layer k
      distBox(model, "dL" + k, "eL" + k, n);
      log.printf("layer %d (%s, %.1f nm): %d elements%n", k, mat[k], d[k], n);
    }
    vbox(model, "e_sup", w, D, D + hsup);
    vbox(model, "e_sub", w, -hsub, 0);
    distBox(model, "d_sup", "e_sup", (int) Math.ceil(hsup / (lamMin / model.param().evaluate("n_sup") / 16)));
    distBox(model, "d_sub", "e_sub", (int) Math.ceil(hsub / (lamMin / model.param().evaluate("n_sub") / 16)));
    model.component("comp1").mesh("mesh1").run();
    try { log.println("Mesh elements: " + model.component("comp1").mesh("mesh1").getNumElem()); } catch (Exception e) { }

    // ------------------------------------------------------------------ study and sweep
    model.study().create("std1");
    model.study("std1").create("freq", "Frequency");
    model.study("std1").feature("freq").set("plist", "c_const/lam0");
    model.sol().create("sol1"); model.sol("sol1").study("std1"); model.sol("sol1").attach("std1"); model.sol("sol1").createAutoSequence("std1");
    model.result().dataset().create("dset_ml", "Solution");
    model.result().dataset("dset_ml").set("solution", "sol1");
    model.result().numerical().create("gev", "EvalGlobal");
    model.result().numerical("gev").set("data", "dset_ml");

    String[] ex = new String[3 + NL]; String[] un = new String[3 + NL];
    ex[0] = "Rr"; ex[1] = "Tt"; ex[2] = "A_Q";
    for (int k = 0; k < NL; k++) ex[3 + k] = "A_L" + k;
    Arrays.fill(un, "1");
    PrintWriter out = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_RTA.txt"));
    StringBuilder hd = new StringBuilder("% lambda_nm R T A_Q 1-R-T");
    for (int k = 0; k < NL; k++) hd.append(" A_L").append(k).append("(").append(mat[k]).append(")");
    out.println(hd);
    for (double L : lams) {
      model.param().set("lam0", L + "[nm]");
      model.sol("sol1").runAll();
      double[] v = evalSafe(model, "gev", ex, un, log);
      StringBuilder sb = new StringBuilder(String.format(Locale.US, "%8.2f %12.6e %12.6e %12.6e %12.6e", L, v[0], v[1], v[2], 1 - v[0] - v[1]));
      for (int k = 0; k < NL; k++) sb.append(String.format(Locale.US, " %12.6e", v[3 + k]));
      out.println(sb); out.flush();
    }
    out.close();

    // ------------------------------------------------------------------ heat-source profile at lam_profile
    model.param().set("lam0", lamProf + "[nm]");
    model.sol("sol1").runAll();
    model.result().dataset().create("cl_z", "CutLine2D");
    model.result().dataset("cl_z").set("data", "dset_ml");
    model.result().dataset("cl_z").set("genpoints", new double[][]{{w / 2, D}, {w / 2, 0}});
    model.result().export().create("ex_q", "Data");
    model.result().export("ex_q").set("data", "cl_z");
    model.result().export("ex_q").set("expr", new String[]{"emw.Qh/I_inc", "abs(emw.Ez)^2+abs(emw.Ex)^2+abs(emw.Ey)^2"});
    model.result().export("ex_q").set("unit", new String[]{"1/m", "V^2/m^2"});
    try { model.result().export("ex_q").set("location", "regulargrid"); model.result().export("ex_q").set("regulargridx1", 2000); } catch (Exception e) { }
    model.result().export("ex_q").set("filename", BASE + "results/" + tag + "_qz.txt");
    model.result().export("ex_q").run();

    model.result().create("pg_E", "PlotGroup2D");
    model.result("pg_E").set("data", "dset_ml");
    model.result("pg_E").create("surf1", "Surface");
    model.result("pg_E").feature("surf1").set("expr", "emw.normE");
    log.println("END");
    log.close();
    return model;
  }

  static void rect(Model m, String tag, String y0, String h) {
    m.component("comp1").geom("geom1").create(tag, "Rectangle");
    m.component("comp1").geom("geom1").feature(tag).set("size", new String[]{"w_cell", h});
    m.component("comp1").geom("geom1").feature(tag).set("pos", new String[]{"0", y0});
  }
  static void box(Model m, String tag, int dim, double x0, double x1, double y0, double y1) {
    m.component("comp1").selection().create(tag, "Box");
    m.component("comp1").selection(tag).set("entitydim", dim);
    m.component("comp1").selection(tag).set("xmin", x0 + ""); m.component("comp1").selection(tag).set("xmax", x1 + "");
    m.component("comp1").selection(tag).set("ymin", y0 + ""); m.component("comp1").selection(tag).set("ymax", y1 + "");
    m.component("comp1").selection(tag).set("condition", "inside");
  }
  static void distBox(Model m, String tag, String sel, int n) {
    m.component("comp1").mesh("mesh1").feature("map1").create(tag, "Distribution");
    m.component("comp1").mesh("mesh1").feature("map1").feature(tag).selection().named(sel);
    m.component("comp1").mesh("mesh1").feature("map1").feature(tag).set("numelem", n);
  }

  static void vbox(Model m, String tag, double w, double y0, double y1) {
    double e = 1e-3;
    box(m, tag + "_l", 1, -e, e, y0 - e, y1 + e);
    box(m, tag + "_r", 1, w - e, w + e, y0 - e, y1 + e);
    m.component("comp1").selection().create(tag, "Union");
    m.component("comp1").selection(tag).set("entitydim", 1);
    m.component("comp1").selection(tag).set("input", new String[]{tag + "_l", tag + "_r"});
  }
  // ---------------------------------------------------------------- shared helpers (kept inline: comsol compile works per file)
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
    if (s.contains(":")) {
      String[] p = s.split(":");
      double a = Double.parseDouble(p[0]), st = Double.parseDouble(p[1]), b = Double.parseDouble(p[2]);
      int n = (int) Math.floor((b - a) / st + 1e-9) + 1;
      double[] o = new double[n]; for (int i = 0; i < n; i++) o[i] = a + i * st; return o;
    }
    String[] p = s.split(",");
    double[] o = new double[p.length]; for (int i = 0; i < p.length; i++) o[i] = Double.parseDouble(p[i].trim()); return o;
  }
  static void material(Model m, String tag, String sel, String eps, String unused) {
    m.component("comp1").material().create(tag, "Common");
    m.component("comp1").material(tag).selection().named(sel);
    String[][] pr = {{"relpermittivity", eps}, {"relpermeability", "1"}, {"electricconductivity", "0"}};
    for (String[] p : pr) m.component("comp1").material(tag).propertyGroup("def").set(p[0], new String[]{p[1], "0", "0", "0", p[1], "0", "0", "0", p[1]});
  }
  static void intop(Model m, String tag, String sel) {
    m.component("comp1").cpl().create(tag, "Integration");
    m.component("comp1").cpl(tag).set("opname", tag);
    m.component("comp1").cpl(tag).selection().named(sel);
  }
  static double[] evalSafe(Model m, String gev, String[] ex, String[] un, PrintWriter log) {
    double[] out = new double[ex.length];
    for (int k = 0; k < ex.length; k++) {
      try {
        m.result().numerical(gev).set("expr", new String[]{ex[k]});
        m.result().numerical(gev).set("unit", new String[]{un[k]});
        out[k] = m.result().numerical(gev).getReal()[0][0];
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
