package prot;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Collects existing stops-v1 comparisons without pooling genomes or reclassifying calls.
 * Rates are reference concordance, not biological precision. @author Sayu
 */
public final class ReferenceConcordanceReport{

	public static void main(String[] args){
		try{run(args);}catch(Throwable t){t.printStackTrace(System.err); System.exit(1);}
	}

	public static void run(String[] args){
		final PreParser pp=new PreParser(args, ReferenceConcordanceReport.class, false);
		String in=null, out=null;
		final HashSet<String> flags=new HashSet<String>();
		for(String arg : pp.args){
			final int eq=arg.indexOf('='); require(eq>0 && eq<arg.length()-1, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			require(flags.add(key), "Duplicate report flag: "+key);
			if(key.equals("in")){in=value;}else if(key.equals("out")){out=value;}
			else{throw new IllegalArgumentException("Unknown report flag: "+key);}
		}
		require(in!=null && out!=null, "ReferenceConcordanceReport requires in=manifest.tsv out=PREFIX");
		final ArrayList<Comparison> comparisons=readManifest(in);
		final LinkedHashMap<String,Arm> arms=new LinkedHashMap<String,Arm>();
		final LinkedHashMap<String,TreeMap<String,String>> references=new LinkedHashMap<String,TreeMap<String,String>>();
		for(Comparison c : comparisons){
			c.summary=readSummary(c.summaryPath);
			final TreeMap<String,String> reference=c.summary.get("reference");
			validateMetrics(reference); validateMetrics(c.summary.get("baseline")); validateMetrics(c.summary.get("candidate"));
			final long r=count(reference, "primary_stops");
			require(count(reference, "reference_primary_stops")==r && count(reference, "matched_stops")==r,
				"Reference inventory must be its own complete intersection: "+c.summaryPath);
			final TreeMap<String,String> old=references.put(c.genome, reference);
			require(old==null || old.equals(reference), "Reference inventory differs between comparisons of "+c.genome);
			addArm(arms, c, "baseline", c.before, c.beforeTime, r);
			addArm(arms, c, "candidate", c.after, c.afterTime, r);
			validateStops(c);
		}
		final String[] outputs={out+".arms.tsv", out+".comparisons.tsv", out+".strata.tsv", out+".changes.tsv"};
		require(Tools.testOutputFiles(false, false, false, outputs), "Report outputs must be new writable files");
		// All input validation precedes opening any output writer.
		writeArms(outputs[0], arms); writeComparisons(outputs[1], comparisons);
		writeStrata(outputs[2], references, arms); writeChanges(outputs[3], comparisons);
		pp.outstream.println("REFERENCE_CONCORDANCE_REPORT_PASS genomes="+references.size()+" arms="+arms.size()
			+" comparisons="+comparisons.size()+" reference_concordance_only=true");
		Shared.closeStream(pp.outstream);
	}

	private static ArrayList<Comparison> readManifest(String path){
		final ArrayList<Comparison> result=new ArrayList<Comparison>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile file=ByteFile.makeByteFile(path, false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(file, lp, MANIFEST);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line); require(lp.terms()==7, "Expected seven manifest columns: "+path);
				final Comparison c=new Comparison(); c.genome=lp.parseString(0); c.before=lp.parseString(1); c.after=lp.parseString(2);
				require(safe(c.genome) && safe(c.before) && safe(c.after) && !c.before.equals(c.after), "Invalid genome/arm identity");
				require(seen.add(c.genome+"\t"+c.before+"\t"+c.after), "Duplicate comparison: "+c.genome+"/"+c.before+"/"+c.after);
				c.summaryPath=lp.parseString(3); c.stopsPath=lp.parseString(4); c.beforeTime=lp.parseString(5); c.afterTime=lp.parseString(6);
				require(Tools.testInputFiles(false, true, c.summaryPath, c.stopsPath, c.beforeTime, c.afterTime), "Unreadable report inputs");
				result.add(c);
			}
		}finally{require(!file.close(), "Manifest read failed: "+path);}
		require(!result.isEmpty(), "Empty report manifest"); return result;
	}

	private static LinkedHashMap<String,TreeMap<String,String>> readSummary(String path){
		final LinkedHashMap<String,TreeMap<String,String>> result=new LinkedHashMap<String,TreeMap<String,String>>();
		for(String arm : new String[]{"reference", "baseline", "candidate", "delta"}){result.put(arm, new TreeMap<String,String>());}
		final ByteFile file=ByteFile.makeByteFile(path, false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(file, lp, SUMMARY);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line); require(lp.terms()==4 && lp.termEquals("stops-v1", 0), "Unknown summary schema: "+path);
				final TreeMap<String,String> metrics=result.get(lp.parseString(1));
				require(metrics!=null, "Unknown summary arm: "+path);
				final String metric=lp.parseString(2), value=lp.parseString(3);
				require(!metric.isEmpty() && !value.isEmpty() && metrics.put(metric, value)==null, "Duplicate/empty summary metric: "+path+"/"+metric);
			}
		}finally{require(!file.close(), "Summary read failed: "+path);}
		for(TreeMap<String,String> metrics : result.values()){require(!metrics.isEmpty(), "Missing summary arm: "+path);}
		return result;
	}

	public static void validateMetrics(TreeMap<String,String> m){
		final long raw=count(m, "raw_cds_rows"), primary=count(m, "primary_rows"), excluded=count(m, "excluded_rows");
		final long c=count(m, "primary_stops"), r=count(m, "reference_primary_stops"), matched=count(m, "matched_stops");
		require(raw==Math.addExact(primary, excluded) && c<=primary, "Raw rows must partition; distinct stops cannot exceed primary rows");
		require(Math.addExact(matched, count(m, "unmatched_primary_stops"))==c
			&& Math.addExact(matched, count(m, "missed_reference_stops"))==r, "Intersection/denominator arithmetic is inconsistent");
		long strata=0;
		for(String key : m.keySet()){if(key.startsWith("rows.")){strata=Math.addExact(strata, count(m, key));}}
		require(strata==raw, "Exclusive row strata must sum to raw CDS rows");
		rate(m, "stop_reference_precision", matched, c); rate(m, "stop_reference_recall", matched, r);
	}

	private static void addArm(LinkedHashMap<String,Arm> arms, Comparison c, String role, String name, String timing, long r){
		final Arm arm=new Arm(); arm.genome=c.genome; arm.name=name; arm.metrics=c.summary.get(role);
		arm.summary=c.summaryPath; arm.role=role; arm.timing=timing; arm.seconds=readTiming(timing);
		require(count(arm.metrics, "reference_primary_stops")==r, "Every arm must use the fixed reference denominator");
		final String key=c.genome+"\t"+name; final Arm old=arms.get(key);
		if(old==null){arms.put(key, arm);}
		else{require(old.metrics.equals(arm.metrics) && old.timing.equals(timing) && old.seconds==arm.seconds,
			"Repeated genome/arm has conflicting evidence: "+key);}
	}

	private static double readTiming(String path){
		final ByteFile file=ByteFile.makeByteFile(path, false); final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[] line=file.nextLine(); require(line!=null, "Empty timing file: "+path); lp.set(line);
			require(lp.terms()==2 && lp.termEquals("elapsed_seconds", 0) && file.nextLine()==null, "Expected one caller elapsed_seconds row: "+path);
			final double value=Double.parseDouble(lp.parseString(1));
			require(Double.isFinite(value) && value>=0, "Invalid caller time: "+path); return value;
		}finally{require(!file.close(), "Timing read failed: "+path);}
	}

	private static void validateStops(Comparison c){
		final long[] delta=new long[DELTA.length], totals=new long[5];
		final HashSet<String> keys=new HashSet<String>();
		final ByteFile file=ByteFile.makeByteFile(c.stopsPath, false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(file, lp, STOPS);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line); require(lp.terms()==14 && lp.termEquals("stops-v1", 0), "Malformed paired-stop row: "+c.stopsPath);
				final String key=lp.parseString(1)+"\t"+lp.parseString(2)+"\t"+lp.parseString(3);
				require(keys.add(key), "Duplicate paired-stop key: "+key);
				for(int i=4; i<=8; i+=2){require(lp.parseLong(i)>=lp.parseLong(i+1) && lp.parseLong(i+1)>=0, "Raw/primary stop multiplicities are inconsistent");}
				final boolean ref=lp.parseLong(5)>0, before=lp.parseLong(7)>0, after=lp.parseLong(9)>0;
				if(ref){totals[0]++;} if(before){totals[1]++; if(ref){totals[3]++;}}
				if(after){totals[2]++; if(ref){totals[4]++;}}
				final boolean changedStart=lp.termEquals("same_stop_changed_starts", 10);
				final String expected=before && after ? (changedStart ? "same_stop_changed_starts" : "unchanged")
					: after ? "gained" : before ? "lost" : "excluded_or_uncalled";
				require(lp.termEquals(expected, 10), "Change label contradicts primary presence: "+key);
				if(before && after){delta[0]++; if(changedStart){delta[5]++;}}
				else if(after){delta[1]++; if(ref){delta[3]++;}}
				else if(before){delta[2]++; if(ref){delta[4]++;}}
				if(changedStart || before!=after){c.changes.add(line.clone());}
			}
		}finally{require(!file.close(), "Paired-stop read failed: "+c.stopsPath);}
		for(int i=0; i<DELTA.length; i++){require(delta[i]==count(c.summary.get("delta"), DELTA[i]), "Paired rows disagree with delta "+DELTA[i]+": "+c.stopsPath);}
		require(totals[0]==count(c.summary.get("reference"), "primary_stops")
			&& totals[1]==count(c.summary.get("baseline"), "primary_stops") && totals[2]==count(c.summary.get("candidate"), "primary_stops")
			&& totals[3]==count(c.summary.get("baseline"), "matched_stops") && totals[4]==count(c.summary.get("candidate"), "matched_stops"),
			"Paired-stop inventories disagree with primary counts/intersections: "+c.stopsPath);
		assert(c.changes.size()==delta[1]+delta[2]+delta[5]) : "Changes output must retain every gain, loss and start-only change";
	}

	private static void writeArms(String path, LinkedHashMap<String,Arm> arms){
		final ByteStreamWriter out=HmmComparisonData.writer(path); final ByteBuilder row=new ByteBuilder();
		try{
			out.println("genome\tarm\tprimary_C\tprimary_R\tmatched_M\tunmatched\tprecision\trecall\tcaller_seconds\traw_cds_rows\texcluded_rows\tsummary\tsummary_role\ttiming");
			for(Arm a : arms.values()){
				row.clear().append(a.genome).tab().append(a.name);
				for(String metric : ARM_METRICS){row.tab().append(a.metrics.get(metric));}
				out.println(row.tab().append(Double.toString(a.seconds)).tab().append(a.metrics.get("raw_cds_rows"))
					.tab().append(a.metrics.get("excluded_rows")).tab().append(a.summary).tab().append(a.role).tab().append(a.timing));
			}
		}finally{HmmComparisonData.close(out);}
	}

	private static void writeComparisons(String path, ArrayList<Comparison> comparisons){
		final ByteStreamWriter out=HmmComparisonData.writer(path); final ByteBuilder row=new ByteBuilder();
		try{
			out.println("genome\tbaseline\tcandidate\tshared_primary_stops\tgained_primary_stops\tlost_primary_stops\tgained_reference_stops\tlost_reference_stops\tsame_stop_changed_starts\tsummary\tstops");
			for(Comparison c : comparisons){
				row.clear().append(c.genome).tab().append(c.before).tab().append(c.after);
				for(String metric : DELTA){row.tab().append(count(c.summary.get("delta"), metric));}
				out.println(row.tab().append(c.summaryPath).tab().append(c.stopsPath));
			}
		}finally{HmmComparisonData.close(out);}
	}

	private static void writeStrata(String path, LinkedHashMap<String,TreeMap<String,String>> refs, LinkedHashMap<String,Arm> arms){
		final ByteStreamWriter out=HmmComparisonData.writer(path); final ByteBuilder row=new ByteBuilder();
		try{
			out.println("genome\trole\tarm\tstratum\trows");
			for(Map.Entry<String,TreeMap<String,String>> e : refs.entrySet()){strata(out, row, e.getKey(), "reference", ".", e.getValue());}
			for(Arm a : arms.values()){strata(out, row, a.genome, "prediction", a.name, a.metrics);}
		}finally{HmmComparisonData.close(out);}
	}
	private static void strata(ByteStreamWriter out, ByteBuilder row, String genome, String role, String arm, TreeMap<String,String> m){
		assert(m!=null) : "Strata are emitted only from already validated summaries";
		for(String key : m.keySet()){if(key.startsWith("rows.")){
			out.println(row.clear().append(genome).tab().append(role).tab().append(arm).tab().append(key.substring(5)).tab().append(count(m, key)));
		}}
	}
	private static void writeChanges(String path, ArrayList<Comparison> comparisons){
		final ByteStreamWriter out=HmmComparisonData.writer(path); final ByteBuilder row=new ByteBuilder();
		try{
			out.println("genome\tbaseline\tcandidate\t"+String.join("\t", STOPS));
			for(Comparison c : comparisons){for(byte[] change : c.changes){
				out.println(row.clear().append(c.genome).tab().append(c.before).tab().append(c.after).tab().append(change));
			}}
		}finally{HmmComparisonData.close(out);}
	}

	private static void header(ByteFile file, LineParser1 lp, String[] names){
		final byte[] line=file.nextLine(); require(line!=null, "Missing TSV header"); lp.set(line);
		require(lp.terms()==names.length, "Unexpected TSV header width");
		for(int i=0; i<names.length; i++){require(lp.termEquals(names[i], i), "Unexpected TSV header at column "+i+": expected "+names[i]);}
	}
	private static long count(Map<String,String> m, String key){
		final String text=m.get(key); require(text!=null, "Missing metric: "+key);
		final long value=Long.parseLong(text); require(value>=0, "Negative count: "+key); return value;
	}
	private static void rate(Map<String,String> m, String key, long numerator, long denominator){
		final String text=m.get(key); require(text!=null, "Missing rate: "+key);
		if(denominator==0){require(text.equals("NA"), "Zero denominator requires NA: "+key);}
		else{final double value=Double.parseDouble(text); require(Double.isFinite(value) && Math.abs(value-numerator/(double)denominator)<=1e-12,
			"Rate disagrees with explicit numerator/denominator: "+key);}
	}
	private static boolean safe(String s){return s.matches("[A-Za-z0-9][A-Za-z0-9_.-]*");}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Arm{String genome, name, summary, role, timing; double seconds; TreeMap<String,String> metrics;}
	private static final class Comparison{
		String genome, before, after, summaryPath, stopsPath, beforeTime, afterTime;
		LinkedHashMap<String,TreeMap<String,String>> summary;
		final ArrayList<byte[]> changes=new ArrayList<byte[]>();
	}
	private static final String[] MANIFEST={"genome", "baseline", "candidate", "summary", "stops", "baseline_timing", "candidate_timing"};
	private static final String[] SUMMARY={"schema", "arm", "metric", "value"};
	private static final String[] STOPS={"schema", "contig", "strand", "terminal", "reference_raw_rows", "reference_primary_rows", "baseline_raw_rows", "baseline_primary_rows", "candidate_raw_rows", "candidate_primary_rows", "change", "reference_rows", "baseline_rows", "candidate_rows"};
	private static final String[] DELTA={"shared_primary_stops", "gained_primary_stops", "lost_primary_stops", "gained_reference_stops", "lost_reference_stops", "same_stop_changed_starts"};
	private static final String[] ARM_METRICS={"primary_stops", "reference_primary_stops", "matched_stops", "unmatched_primary_stops", "stop_reference_precision", "stop_reference_recall"};
	private ReferenceConcordanceReport(){}
}
