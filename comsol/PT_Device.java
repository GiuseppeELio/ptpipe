/*
 * PT_Device.java  -  ptpipe COMSOL template, level L4 (device scale, mm, s)
 *
 * 2D axisymmetric substrate (radius Rd, thickness d_sub) coated by a thin absorbing stack
 * (thickness d_st, explicit domain with a mapped mesh) illuminated by a Gaussian CW beam
 * switched on and off. The heat source in the stack is
 *     q(r, t) = A_abs * I(r) * rect(t) / d_st,     I(r) = 2 P_inc/(pi w0^2) exp(-2 r^2/w0^2)
 * with A_abs taken from the optical level (PT_Multilayer / TMM, or n_s*sigma_abs for a
 * particle monolayer). Convective + radiative losses on all external boundaries and a closed
 * energy budget  P_abs = dE/dt + Q_conv + Q_rad  (as in Lio et al., Adv. Funct. Mater. 2026).
 *
 * Configuration: comsol/device.cfg     Output: results/<tag>_timeseries.txt, <tag>_T_r.txt
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class PT_Device {

  static String BASE = "/path/to/ptpipe/";   // set automatically by run_validation.sh

  public static Model run() {
    try { return runImpl(); } catch (Exception e) { throw new RuntimeException(e); }
  }

  static Model runImpl() throws Exception {
    try { String r = System.getProperty("ptpipe.root"); if (r != null) BASE = r.endsWith("/") ? r : r + "/"; } catch (Exception e) { }
    Map<String, String> cfg = readCfg(BASE + "comsol/device.cfg");
    String tag = g(cfg, "tag", "device");
    PrintWriter log = new PrintWriter(new FileWriter(BASE + "results/" + tag + "_log.txt"), true);
    log.println("cfg: " + cfg);

    Model model = ModelUtil.create("Model");
    model.modelPath(BASE + "comsol");
    model.label("PT_Device_" + tag + ".mph");

    // ------------------------------------------------------------------ parameters
    model.param().set("A_abs", g(cfg, "A_abs", "0.884"), "Absorbed fraction at the laser wavelength (optical level)");
    model.param().set("P_inc", g(cfg, "P_inc_W", "0.2") + "[W]", "Laser power");
    model.param().set("w0", g(cfg, "w0_mm", "0.5") + "[mm]", "1/e^2 beam radius");
    model.param().set("I_peak", "2*P_inc/(pi*w0^2)");
    model.param().set("t_on", g(cfg, "t_on", "2") + "[s]");
    model.param().set("t_off", g(cfg, "t_off", "122") + "[s]");
    model.param().set("t_stop", g(cfg, "t_stop", "242") + "[s]");
    model.param().set("dt_out", g(cfg, "dt_out", "0.1") + "[s]");
    model.param().set("dt_s", g(cfg, "dt_smooth", "0.05") + "[s]", "Smoothing of the on/off edges");
    model.param().set("Rd", g(cfg, "Rd_mm", "12.5") + "[mm]");
    model.param().set("d_sub", g(cfg, "d_sub_mm", "1.0") + "[mm]");
    model.param().set("d_st", g(cfg, "d_st_nm", "230") + "[nm]");
    model.param().set("z_top", "d_sub + d_st");
    model.param().set("k_sub", g(cfg, "k_sub", "1.38") + "[W/(m*K)]");
    model.param().set("rho_sub", g(cfg, "rho_sub", "2203") + "[kg/m^3]");
    model.param().set("cp_sub", g(cfg, "cp_sub", "703") + "[J/(kg*K)]");
    model.param().set("k_st", g(cfg, "k_st", "168.6") + "[W/(m*K)]", "Stack conductivity (thickness-weighted, in-plane)");
    model.param().set("rho_st", g(cfg, "rho_st", "15000") + "[kg/m^3]");
    model.param().set("cp_st", g(cfg, "cp_st", "200") + "[J/(kg*K)]");
    model.param().set("h_top", g(cfg, "h_top", "10") + "[W/(m^2*K)]");
    model.param().set("h_bot", g(cfg, "h_bot", "10") + "[W/(m^2*K)]");
    model.param().set("h_edge", g(cfg, "h_edge", "10") + "[W/(m^2*K)]");
    model.param().set("emis", g(cfg, "emis", "0.9"));
    model.param().set("T_amb", "293.15[K]");

    model.func().create("rect1", "Rectangle");
    model.func("rect1").set("lower", "t_on");
    model.func("rect1").set("upper", "t_off");
    model.func("rect1").set("smooth", "dt_s");

    // ------------------------------------------------------------------ geometry
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 2);
    model.component("comp1").geom("geom1").axisymmetric(true);
    model.component("comp1").geom("geom1").create("r_sub", "Rectangle");
    model.component("comp1").geom("geom1").feature("r_sub").set("size", new String[]{"Rd", "d_sub"});
    model.component("comp1").geom("geom1").feature("r_sub").set("pos", new String[]{"0", "0"});
    model.component("comp1").geom("geom1").create("r_st", "Rectangle");
    model.component("comp1").geom("geom1").feature("r_st").set("size", new String[]{"Rd", "d_st"});
    model.component("comp1").geom("geom1").feature("r_st").set("pos", new String[]{"0", "d_sub"});
    model.component("comp1").geom("geom1").run("fin");

    // ------------------------------------------------------------------ selections
    String e = "0.1*d_st";      // 23 nm for the default stack: well above the geometric tolerance
    box(model, "sel_st", 2, "-" + e, "Rd+" + e, "d_sub-" + e, "z_top+" + e, "Stack");
    box(model, "sel_sub", 2, "-" + e, "Rd+" + e, "-" + e, "d_sub+" + e, "Substrate");
    box(model, "sel_top", 1, "-" + e, "Rd+" + e, "z_top-" + e, "z_top+" + e, "Top surface");
    box(model, "sel_bot", 1, "-" + e, "Rd+" + e, "-" + e, e, "Bottom surface");
    box(model, "sel_edge", 1, "Rd-" + e, "Rd+" + e, "-" + e, "z_top+" + e, "Outer edge");
    box(model, "sel_ax_st", 1, "-" + e, e, "d_sub-" + e, "z_top+" + e, "Axis in stack");
    box(model, "sel_ax_sub", 1, "-" + e, e, "-" + e, "d_sub+" + e, "Axis in substrate");
    box(model, "sel_int", 1, "-" + e, "Rd+" + e, "d_sub-" + e, "d_sub+" + e, "Stack/substrate interface");
    model.component("comp1").selection().create("sel_horiz", "Union");
    model.component("comp1").selection("sel_horiz").set("entitydim", 1);
    model.component("comp1").selection("sel_horiz").set("input", new String[]{"sel_top", "sel_int", "sel_bot"});
    model.component("comp1").selection().create("sel_ext", "Union");
    model.component("comp1").selection("sel_ext").set("entitydim", 1);
    model.component("comp1").selection("sel_ext").set("input", new String[]{"sel_top", "sel_bot", "sel_edge"});
    box(model, "pt_top_s", 0, "-" + e, e, "z_top-" + e, "z_top+" + e, "Top centre");
    box(model, "pt_bot_s", 0, "-" + e, e, "-" + e, e, "Bottom centre");

    try {
      String[] ss = {"sel_st", "sel_sub", "sel_top", "sel_bot", "sel_edge", "sel_ax_st", "sel_ax_sub", "sel_int", "pt_top_s", "pt_bot_s"};
      int[] dd = {2, 2, 1, 1, 1, 1, 1, 1, 0, 0};
      for (int q = 0; q < ss.length; q++) {
        int[] en = model.component("comp1").selection(ss[q]).entities(dd[q]);
        log.println("Selection " + ss[q] + ": " + Arrays.toString(en));
        if (en.length == 0) log.println("WARN empty selection " + ss[q]);
      }
    } catch (Exception ex) { log.println("WARN diag: " + ex.getMessage()); }

    // ------------------------------------------------------------------ variables
    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("I_r", "I_peak*exp(-2*(r/w0)^2)");
    model.component("comp1").variable("var1").set("Q_las", "A_abs*I_r*rect1(t)/d_st", "Volumetric heat source in the stack");
    model.component("comp1").variable("var1").set("Q_stored", "intop_st(rho_st*cp_st*(T-T_amb)) + intop_sub(rho_sub*cp_sub*(T-T_amb))");
    model.component("comp1").variable("var1").set("Q_abs", "intop_st(Q_las)");
    model.component("comp1").variable("var1").set("Q_conv", "intop_top(h_top*(T-T_amb)) + intop_bot(h_bot*(T-T_amb)) + intop_edge(h_edge*(T-T_amb))");
    model.component("comp1").variable("var1").set("Q_rad", "intop_ext(emis*sigma_const*(T^4-T_amb^4))");
    model.component("comp1").variable("var1").set("T_top_c", "pt_top(T)");
    model.component("comp1").variable("var1").set("T_bot_c", "pt_bot(T)");
    intop(model, "intop_st", 2, "sel_st", true);
    intop(model, "intop_sub", 2, "sel_sub", true);
    intop(model, "intop_top", 1, "sel_top", true);
    intop(model, "intop_bot", 1, "sel_bot", true);
    intop(model, "intop_edge", 1, "sel_edge", true);
    intop(model, "intop_ext", 1, "sel_ext", true);
    intop(model, "pt_top", 0, "pt_top_s", false);
    intop(model, "pt_bot", 0, "pt_bot_s", false);

    // ------------------------------------------------------------------ materials
    mat(model, "mat_st", "sel_st", "k_st", "rho_st", "cp_st", "Absorbing stack");
    mat(model, "mat_sub", "sel_sub", "k_sub", "rho_sub", "cp_sub", "Substrate");

    // ------------------------------------------------------------------ physics
    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "T_amb");
    model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
    model.component("comp1").physics("ht").feature("hs1").selection().named("sel_st");
    model.component("comp1").physics("ht").feature("hs1").set("Q0", "Q_las");
    conv(model, "hf_top", "sel_top", "h_top");
    conv(model, "hf_bot", "sel_bot", "h_bot");
    conv(model, "hf_edge", "sel_edge", "h_edge");
    model.component("comp1").physics("ht").create("sar1", "SurfaceToAmbientRadiation", 1);
    model.component("comp1").physics("ht").feature("sar1").selection().named("sel_ext");
    model.component("comp1").physics("ht").feature("sar1").set("epsilon_rad_mat", "userdef");
    model.component("comp1").physics("ht").feature("sar1").set("epsilon_rad", "emis");
    try { model.component("comp1").physics("ht").feature("sar1").set("Tamb", "T_amb"); } catch (Exception ex) { }

    // ------------------------------------------------------------------ mapped mesh, uniform distributions
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").create("map1", "Map");
    dist(model, "dis_st", "sel_ax_st", Integer.parseInt(g(cfg, "n_st", "2")));
    dist(model, "dis_sub", "sel_ax_sub", Integer.parseInt(g(cfg, "n_sub", "40")));
    dist(model, "dis_r", "sel_horiz", Integer.parseInt(g(cfg, "n_r", "250")));
    model.component("comp1").mesh("mesh1").run();
    try { log.println("Mesh elements: " + model.component("comp1").mesh("mesh1").getNumElem()); } catch (Exception ex) { }

    // ------------------------------------------------------------------ study
    model.study().create("std1");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("time").set("tlist", "range(0,dt_out,t_stop)");
    try { model.study("std1").feature("time").set("usertol", true); model.study("std1").feature("time").set("rtol", "1e-4"); } catch (Exception ex) { }
    model.sol().create("sol1"); model.sol("sol1").study("std1"); model.sol("sol1").attach("std1"); model.sol("sol1").createAutoSequence("std1");
    try {
      model.sol("sol1").feature("t1").set("maxstepconstraintbdf", "const");
      model.sol("sol1").feature("t1").set("maxstepbdf", "dt_out");
    } catch (Exception ex) { log.println("WARN max step: " + ex.getMessage()); }
    model.sol("sol1").runAll();

    // ------------------------------------------------------------------ results
    model.result().dataset().create("dset_m", "Solution");
    model.result().dataset("dset_m").set("solution", "sol1");
    model.result().table().create("tbl_ts", "Table");
    model.result().numerical().create("gev_ts", "EvalGlobal");
    model.result().numerical("gev_ts").set("data", "dset_m");
    model.result().numerical("gev_ts").set("expr", new String[]{"T_top_c-T_amb", "T_bot_c-T_amb", "Q_abs", "d(Q_stored,t)", "Q_conv", "Q_rad",
        "Q_abs - d(Q_stored,t) - Q_conv - Q_rad"});
    model.result().numerical("gev_ts").set("unit", new String[]{"K", "K", "W", "W", "W", "W", "W"});
    model.result().numerical("gev_ts").set("table", "tbl_ts");
    model.result().numerical("gev_ts").setResult();
    model.result().table("tbl_ts").save(BASE + "results/" + tag + "_timeseries.txt");

    double dto = model.param().evaluate("dt_out"), toff = model.param().evaluate("t_off");
    int lev = (int) Math.round(toff / dto) + 1;
    model.result().dataset().create("cln_r", "CutLine2D");
    model.result().dataset("cln_r").set("data", "dset_m");
    double zt = model.param().evaluate("z_top"), rd = model.param().evaluate("Rd");
    model.result().dataset("cln_r").set("genpoints", new double[][]{{0, zt}, {rd, zt}});
    model.result().create("pg_r", "PlotGroup1D");
    model.result("pg_r").set("data", "cln_r");
    model.result("pg_r").setIndex("looplevelinput", "manual", 0);
    model.result("pg_r").setIndex("looplevel", new int[]{lev}, 0);
    model.result("pg_r").create("lg1", "LineGraph");
    model.result("pg_r").feature("lg1").set("expr", "T-T_amb");
    model.result("pg_r").feature("lg1").set("unit", "K");
    model.result("pg_r").feature("lg1").set("xdata", "expr");
    model.result("pg_r").feature("lg1").set("xdataexpr", "r");
    model.result("pg_r").feature("lg1").set("xdataunit", "mm");
    model.result().export().create("exp_r", "Plot");
    model.result().export("exp_r").set("plotgroup", "pg_r");
    model.result().export("exp_r").set("plot", "lg1");
    model.result().export("exp_r").set("filename", BASE + "results/" + tag + "_T_r.txt");
    model.result().export("exp_r").run();

    model.result().create("pg_T", "PlotGroup2D");
    model.result("pg_T").set("data", "dset_m");
    model.result("pg_T").create("surf1", "Surface");
    model.result("pg_T").feature("surf1").set("expr", "T-T_amb");
    model.result("pg_T").feature("surf1").set("colortable", "HeatCameraLight");
    // the full transient solution is large (GB for thousands of output times): keep it only on request
    if (!g(cfg, "keep_solution", "false").equalsIgnoreCase("true")) {
      try { model.sol("sol1").clearSolutionData(); log.println("solution data cleared before saving (keep_solution = true to keep it)"); }
      catch (Exception ex) { log.println("WARN clearSolutionData: " + ex.getMessage()); }
    }
    log.println("END");
    log.close();
    return model;
  }

  // ---------------------------------------------------------------- helpers
  static Map<String, String> readCfg(String f) throws IOException {
    Map<String, String> c = new LinkedHashMap<>();
    BufferedReader br = new BufferedReader(new FileReader(f));
    String s;
    while ((s = br.readLine()) != null) {
      int h = s.indexOf('#'); if (h >= 0) s = s.substring(0, h);
      int q = s.indexOf('='); if (q < 0) continue;
      c.put(s.substring(0, q).trim(), s.substring(q + 1).trim());
    }
    br.close();
    return c;
  }
  static String g(Map<String, String> c, String k, String d) { return c.containsKey(k) ? c.get(k) : d; }
  static void box(Model m, String tag, int dim, String x0, String x1, String y0, String y1, String label) {
    m.component("comp1").selection().create(tag, "Box");
    m.component("comp1").selection(tag).label(label);
    m.component("comp1").selection(tag).set("entitydim", dim);
    m.component("comp1").selection(tag).set("xmin", x0); m.component("comp1").selection(tag).set("xmax", x1);
    m.component("comp1").selection(tag).set("ymin", y0); m.component("comp1").selection(tag).set("ymax", y1);
    m.component("comp1").selection(tag).set("condition", "inside");
  }
  static void intop(Model m, String tag, int dim, String sel, boolean axisym) {
    m.component("comp1").cpl().create(tag, "Integration");
    m.component("comp1").cpl(tag).set("opname", tag);
    m.component("comp1").cpl(tag).set("axisym", axisym);
    m.component("comp1").cpl(tag).selection().named(sel);
  }
  static void mat(Model m, String tag, String sel, String k, String rho, String cp, String label) {
    m.component("comp1").material().create(tag, "Common");
    m.component("comp1").material(tag).label(label);
    m.component("comp1").material(tag).selection().named(sel);
    m.component("comp1").material(tag).propertyGroup("def").set("thermalconductivity", new String[]{k, "0", "0", "0", k, "0", "0", "0", k});
    m.component("comp1").material(tag).propertyGroup("def").set("density", rho);
    m.component("comp1").material(tag).propertyGroup("def").set("heatcapacity", cp);
  }
  static void conv(Model m, String tag, String sel, String h) {
    m.component("comp1").physics("ht").create(tag, "HeatFluxBoundary", 1);
    m.component("comp1").physics("ht").feature(tag).selection().named(sel);
    m.component("comp1").physics("ht").feature(tag).set("HeatFluxType", "ConvectiveHeatFlux");
    m.component("comp1").physics("ht").feature(tag).set("h", h);
    m.component("comp1").physics("ht").feature(tag).set("Text", "T_amb");
  }
  static void dist(Model m, String tag, String sel, int n) {
    m.component("comp1").mesh("mesh1").feature("map1").create(tag, "Distribution");
    m.component("comp1").mesh("mesh1").feature("map1").feature(tag).selection().named(sel);
    m.component("comp1").mesh("mesh1").feature("map1").feature(tag).set("numelem", n);
  }

  public static void main(String[] args) { run(); }
}
