package bin;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;

/** Focused report-import regressions; no ProkCC runtime or model dependency.
 * @author Collei
 */
public final class GradeBinsProkCCTest {

	/** Prepares fixtures or verifies the completed GradeBins CLI reports. */
	public static void main(String[] args){
		String mode=null, dir=null, real=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("mode")){mode=value;}
			else if(key.equals("dir")){dir=value;}
			else if(key.equals("real")){real=value;}
			else{throw new IllegalArgumentException("Unknown test flag: "+key);}
		}
		require(dir!=null, "A private fixture directory is required");
		if("prepare".equals(mode)){prepare(dir, real);}
		else if("verify".equals(mode)){verify(dir);}
		else if("rna".equals(mode)){verifyRna(dir);}
		else{throw new IllegalArgumentException("mode must be prepare or verify");}
		System.out.println("GRADEBINS_PROKCC_"+mode.toUpperCase(java.util.Locale.ROOT)+"_PASS");
	}

	/** Checks real producer output and malformed-input behavior before creating CLI fixtures. */
	private static void prepare(String dir, String real){
		require(new File(dir).mkdir(), "Fixture directory must be new: "+dir);
		require(real!=null, "The native ProkCC report is required for round-trip validation");
		final Map<String, ?> actual=GradeBins.loadProkCC(real);
		final HashMap<String, String> expected=new HashMap<String, String>();
		int rows=0;
		//Independent fixed-column oracle for the inspected v1.2.1 native report.
		for(byte[] bytes : ByteFile.toLines(real)){
			final String line=new String(bytes, java.nio.charset.StandardCharsets.UTF_8);
			if(line.startsWith("#") || line.isEmpty()){continue;}
			final String[] fields=line.split("\t", -1);
			require(fields.length>=4, "Native report row is truncated");
			expected.put(ReadWrite.stripToCore(fields[1]), "-1, "+Float.parseFloat(fields[2])+", "+Float.parseFloat(fields[3]));
			rows++;
		}
		require(rows>0 && actual.size()==expected.size(), "Real report row/key coverage differs");
		for(Map.Entry<String, String> entry : expected.entrySet()){
			require(entry.getValue().equals(String.valueOf(actual.get(entry.getKey()))), "Real fractions/name changed: "+entry.getKey());
		}
		System.out.println("REAL_PROKCC_ROWS\t"+rows);
		final String header="#columns\tgene_contamination\tbin_id\tgene_completeness\textra\n";
		write(dir+"/prokcc.tsv", header+"0.2\t/old/bin.fa\t0.8\tNA\n3e-2\tC:\\bins\\bin.fa.gz\t0.95\tNA\n0\textra.fa\t1\tNA\n");
		final Map<String, ?> named=GradeBins.loadProkCC(dir);
		require(named.size()==2 && "-1, 0.95, 0.03".equals(String.valueOf(named.get("bin"))),
			"Directory, reordered columns, Windows/gzip name or last-duplicate semantics changed");
		final String[] bad={"", "bin\t.5\t.1\n", "#columns\tbin_id\tgene_completeness\n",
			header+"0.1\tbin\t0.9\n", header+"0.1\t\t0.9\tNA\n",
			header+"NaN\tbin\t0.9\tNA\n", header+"Infinity\tbin\t0.9\tNA\n",
			header+"0.1\tbin\t95\tNA\n", header+"0.1\tbin\t-1\tNA\n",
			header+"0.1\tbin\t0.9oops\tNA\n", header+header,
			"#columns\tbin_id\tgene_completeness\tgene_contamination\tbin_id\n"};
		for(int i=0; i<bad.length; i++){
			final String path=dir+"/bad"+i+".tsv";
			write(path, bad[i]);
			boolean rejected=false;
			try{GradeBins.loadProkCC(path);}catch(IllegalArgumentException e){rejected=true;}
			require(rejected, "Malformed ProkCC report accepted: "+i);
		}
		write(dir+"/bin.fa", ">unlabeled\nACGTACGTACGTACGTACGTACGTACGTACGT\n");
		write(dir+"/checkm.tsv", "Name\tCompleteness\tContamination\n/old/bin.fa\t80\t20\nbin.fa.gz\t95\t3\nextra.fa\t100\t0\n");
		write(dir+"/euk-high.tsv", "bin\tcompleteness\tcontamination\nbin.fa\t98\t2\n");
		write(dir+"/euk-low.tsv", "bin\tcompleteness\tcontamination\nbin.fa\t94\t2\n");
		write(dir+"/euk-tie.tsv", "bin\tcompleteness\tcontamination\nbin.fa\t95\t2\n");
		write(dir+"/euk-zero.tsv", "bin\tcompleteness\tcontamination\nbin.fa\t0\t2\n");
		write(dir+"/missing.tsv", header+"0\tother.fa\t1\tNA\n");
		write(dir+"/missing-checkm.tsv", "Name\tCompleteness\tContamination\nother.fa\t100\t0\n");
		final String rnaHeader="#columns\tbin_id\tgene_completeness\tgene_contamination\tr16\tr23\tr5\ttrna\tr18\n";
		write(dir+"/rna.tsv", rnaHeader+"bin.fa\t0.95\t0.03\t2\t3\t4\t20\tNA\n");
		write(dir+"/rna-partial.tsv", rnaHeader.trim()+"\tcds\n"+"bin.fa\t0.95\t0.03\t2\t3\t4\tNA\tNA\t17\n");
		write(dir+"/rna-old.tsv", "#columns\tbin_id\tgene_completeness\tgene_contamination\tr16\tr23\tr5\ttrna\n"
			+"bin.fa\t0.95\t0.03\t2\t3\t4\t20\n");
		for(String count : new String[]{"-1", "1.5", "2147483648", "oops"}){
			write(dir+"/bad-count-"+count+".tsv", rnaHeader+"bin.fa\t0.95\t0.03\t2\t3\t4\t"+count+"\tNA\n");
			boolean rejected=false;
			try{GradeBins.loadProkCC(dir+"/bad-count-"+count+".tsv");}catch(IllegalArgumentException e){rejected=true;}
			require(rejected, "Malformed optional RNA count accepted: "+count);
		}
	}

	/** Checks the actual CLI result, including existing CheckM behavior and absent-bin fallback. */
	private static void verify(String dir){
		checkScore(dir+"/prokcc-out.tsv", "ProkCC", .95f, .03f);
		checkScore(dir+"/checkm-out.tsv", "CheckM2", .95f, .03f);
		checkScore(dir+"/euk-high-out.tsv", "EukCC", .98f, .02f);
		checkScore(dir+"/euk-low-out.tsv", "ProkCC", .95f, .03f);
		checkScore(dir+"/euk-tie-out.tsv", "ProkCC", .95f, .03f);
		checkScore(dir+"/euk-zero-out.tsv", "EukCC", 0, .02f);
		compareGradeCore(dir+"/prokcc-out.tsv", dir+"/checkm-out.tsv", true);
		compareGradeCore(dir+"/missing-out.tsv", dir+"/missing-checkm-out.tsv", false);
		require(report(dir+"/checkm-out.tsv").equals(report(dir+"/baseline-checkm-out.tsv")),
			"CheckM-only output changed from the pre-edit baseline");
	}

	/** ProkCC adds optional count fields, but the existing report fields must remain identical. */
	private static void compareGradeCore(String actual, String expected, boolean prokSource){
		final Map<String, String> a=row(actual, false), e=row(expected, false);
		for(Map.Entry<String, String> entry : e.entrySet()){
			final String value=prokSource && entry.getKey().equals("Source") ? "ProkCC" : entry.getValue();
			require(value.equals(a.get(entry.getKey())), "Imported report changed "+entry.getKey());
		}
	}

	/** Verifies native before/after predictions, same-bin caller counts, and import/unknown behavior. */
	private static void verifyRna(String dir){
		final Map<String, String> before=row(dir+"/standalone-before.tsv", true);
		final Map<String, String> after=row(dir+"/standalone-after.tsv", true);
		for(Map.Entry<String, String> entry : before.entrySet()){
			if(entry.getKey().equals("bin_worker_wall_seconds")){continue;}
			require(entry.getValue().equals(after.get(entry.getKey())), "Standalone ProkCC changed "+entry.getKey());
		}
		require("NA".equals(after.get("r18")), "Disabled ProkCC 18S must be NA");
		final Map<String, String> called=row(dir+"/real-called.tsv", false), imported=row(dir+"/real-imported.tsv", false);
		final String[] nativeNames={"r16", "r23", "r5", "trna"}, gradeNames={"16S", "23S", "5S", "tRNA"};
		for(int i=0; i<nativeNames.length; i++){
			require(after.get(nativeNames[i]).equals(imported.get(gradeNames[i])), "Report import changed "+nativeNames[i]);
			require(called.get(gradeNames[i]).equals(imported.get(gradeNames[i])), "Same-bin gene count differs for "+nativeNames[i]);
		}
		require("NA".equals(imported.get("18S")), "Unknown 18S became a false zero");
		for(String file : new String[]{"rna-out.tsv", "rna-old-out.tsv"}){
			final Map<String, String> values=row(dir+"/"+file, false);
			require("2".equals(values.get("16S")) && "3".equals(values.get("23S")) && "4".equals(values.get("5S"))
				&& "20".equals(values.get("tRNA")) && "NA".equals(values.get("18S")), "Known/unknown imported count changed: "+file);
		}
		final Map<String, String> partial=row(dir+"/rna-partial-out.tsv", false);
		require("2".equals(partial.get("16S")) && "0".equals(partial.get("tRNA")), "Partial report must preserve known values and call for missing MIMAG counts");
		require("17".equals(partial.get("CDS")) && "NA".equals(partial.get("CDSLen")), "Imported CDS counts cannot use another annotation's length");
	}

	/** Independent header-keyed reader for exactly one completed CLI output row. */
	private static Map<String, String> row(String path, boolean prok){
		String[] header=null, data=null;
		for(byte[] bytes : ByteFile.toLines(path)){
			final String line=new String(bytes, java.nio.charset.StandardCharsets.UTF_8);
			if(line.startsWith(prok ? "#columns\t" : "#Num\t")){header=line.split("\t", -1);}
			else if(!line.startsWith("#") && !line.isEmpty()){
				require(data==null, "Expected one row: "+path); data=line.split("\t", -1);
			}
		}
		final int offset=prok ? 1 : 0;
		require(header!=null && data!=null && header.length==data.length+offset, "Report schema mismatch: "+path);
		final Map<String, String> map=new HashMap<String, String>();
		for(int i=0; i<data.length; i++){map.put(header[i+offset], data[i]);}
		return map;
	}

	/** Locates public output columns independently of their optional neighbors. */
	private static void checkScore(String path, String source, float comp, float contam){
		final ArrayList<byte[]> lines=ByteFile.toLines(path);
		require(lines.size()==2, "Expected one bin plus header: "+path);
		final LineParser1 lp=new LineParser1('\t');
		lp.set(lines.get(0));
		int c=-1, t=-1, s=-1;
		for(int i=0; i<lp.terms(); i++){
			if(lp.termEquals("Completeness", i)){c=i;}
			if(lp.termEquals("Contam", i)){t=i;}
			if(lp.termEquals("Source", i)){s=i;}
		}
		require(c>=0 && t>=0 && s>=0, "Missing score/source columns: "+path);
		lp.set(lines.get(1));
		require(lp.parseFloat(c)==comp && lp.parseFloat(t)==contam && lp.termEquals(source, s), "Wrong imported score/source: "+path);
	}

	private static String report(String path){
		final StringBuilder text=new StringBuilder();
		for(byte[] line : ByteFile.toLines(path)){text.append(new String(line, java.nio.charset.StandardCharsets.UTF_8)).append('\n');}
		return text.toString();
	}

	private static void write(String path, String text){
		final ByteStreamWriter out=new ByteStreamWriter(path, false, false, false);
		out.start(); out.print(text);
		require(!out.poisonAndWait(), "Fixture write failed: "+path);
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
