/*
 * PT_Figures.java  -  ptpipe: export COMSOL field maps and meshes of solved cases (no solving).
 *
 * Reads the case names listed in comsol/figures.cfg (one per line) and, for each
 * results/<case>_Model.mph, writes to results/:
 *   <case>_mesh.mphtxt                 mesh (all templates)
 *   PT_Particle:   <case>_E2_zx.txt, <case>_E2_xy.txt   |E|^2/E0^2 on the planes y = 0 and z = 0
 *                  <case>_T_zx.txt                       temperature rise on the plane y = 0
 *   PT_Multilayer: <case>_field2d.txt                    x y |E| (V/m) q/I_inc (1/m)
 * Plots: python validation/plot_fields.py <case> ...
 * Run:   bash run_validation.sh figures      (after the cases have been solved)
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

public class PT_Figures {

  static String BASE = "/path/to/ptpipe/";   // set automatically by run_validation.sh

  public static Model run() {
    try { String r = System.getProperty("ptpipe.root"); if (r != null) BASE = r.endsWith("/") ? r : r + "/"; } catch (Exception e) { }
    List<String> cases = new ArrayList<>();
    try {
      BufferedReader br = new BufferedReader(new FileReader(BASE + "comsol/figures.cfg"));
      String s;
      while ((s = br.readLine()) != null) {
        int h = s.indexOf('#'); if (h >= 0) s = s.substring(0, h);
        s = s.trim(); if (!s.isEmpty()) cases.add(s);
      }
      br.close();
    } catch (Exception e) { System.out.println("ERROR reading comsol/figures.cfg: " + e); return null; }

    PrintWriter log = null;
    try { log = new PrintWriter(new FileWriter(BASE + "results/figures_log.txt", true), true); } catch (Exception e) { }
    for (String c : cases) {
      String mph = BASE + "results/" + c + "_Model.mph";
      if (!new File(mph).exists()) mph = BASE + "results/" + c + ".mph";
      Model m = null;
      try {
        m = ModelUtil.load("F" + Math.abs(c.hashCode()), mph);
        try { m.component("comp1").mesh("mesh1").export(BASE + "results/" + c + "_mesh.mphtxt"); }
        catch (Exception e) { msg(log, c + ": WARN mesh export " + e.getMessage()); }
        List<String> ds = Arrays.asList(m.result().dataset().tags());
        if (ds.contains("dset_em")) {
          plane(m, log, c, "dset_em", "zx", "E2n/E0^2", "1", "_E2_zx.txt");
          plane(m, log, c, "dset_em", "xy", "E2n/E0^2", "1", "_E2_xy.txt");
          if (ds.contains("dset_ht")) plane(m, log, c, "dset_ht", "zx", "T-T0", "K", "_T_zx.txt");
        } else if (ds.contains("dset_ml")) {
          m.result().export().create("fx_ml", "Data");
          m.result().export("fx_ml").set("data", "dset_ml");
          m.result().export("fx_ml").set("expr", new String[]{"emw.normE", "emw.Qh/I_inc"});
          m.result().export("fx_ml").set("unit", new String[]{"V/m", "1/m"});
          m.result().export("fx_ml").set("filename", BASE + "results/" + c + "_field2d.txt");
          m.result().export("fx_ml").run();
          double Iinc = m.param().evaluate("P_in") / (m.param().evaluate("w_cell") * 1.0);
          msg(log, c + ": field2d written; I_inc = " + Iinc + " W/m^2 (per 1 m out of plane)");
        } else {
          msg(log, c + ": no field datasets (device maps are written by PT_Device itself)");
        }
        msg(log, c + ": OK");
      } catch (Exception e) { msg(log, c + ": ERROR " + e); }
      try { if (m != null) ModelUtil.remove(m.tag()); } catch (Exception e) { }
    }
    if (log != null) log.close();
    return null;   // nothing to save
  }

  static void plane(Model m, PrintWriter log, String c, String data, String plane, String expr, String unit, String suffix) {
    String tag = "cp_" + data + "_" + plane;
    try {
      m.result().dataset().create(tag, "CutPlane");
      m.result().dataset(tag).set("data", data);
      m.result().dataset(tag).set("quickplane", plane);
      m.result().export().create("fx_" + tag, "Data");
      m.result().export("fx_" + tag).set("data", tag);
      m.result().export("fx_" + tag).set("expr", new String[]{expr});
      m.result().export("fx_" + tag).set("unit", new String[]{unit});
      m.result().export("fx_" + tag).set("filename", BASE + "results/" + c + suffix);
      m.result().export("fx_" + tag).run();
    } catch (Exception e) { msg(log, c + ": WARN " + suffix + " " + e.getMessage()); }
  }

  static void msg(PrintWriter log, String s) { System.out.println(s); if (log != null) log.println(s); }

  public static void main(String[] a) { run(); }
}
