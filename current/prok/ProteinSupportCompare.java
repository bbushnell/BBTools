package prok;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteStreamWriter;
import gff.GffLine;
import structures.ByteBuilder;

/** Compares CDS coordinates to external GFFs and records every changed call.
 * Exact/stop concordance is reference agreement, not proof of biological truth.
 * @author Keqing
 */
public final class ProteinSupportCompare {

	public static void main(String[] args){
		for(String arg : args){
			if(arg.regionMatches(true, 0, "mode=", 0, 5)){ProteinStopEvaluation.run(args); return;}
		}
		String baseline=null, candidate=null, reference=null, out=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected baseline=GFF candidate=GFF ref=GFF out=TSV");}
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("baseline")){baseline=value;}else if(key.equals("candidate")){candidate=value;}
			else if(key.equals("ref")){reference=value;}else if(key.equals("out")){out=value;}
			else{throw new IllegalArgumentException("Unknown comparison option: "+key);}
		}
		if(baseline==null || candidate==null || reference==null || out==null){
			throw new IllegalArgumentException("Required: baseline=GFF candidate=GFF ref=GFF out=TSV");
		}
		if(!shared.Tools.testOutputFiles(false, false, false, out)){throw new IllegalArgumentException("Output already exists or is unwritable: "+out);}
		final ArrayList<GffLine> before=GffLine.loadGffFile(baseline, "CDS", false);
		final ArrayList<GffLine> after=GffLine.loadGffFile(candidate, "CDS", false);
		final ArrayList<GffLine> truth=GffLine.loadGffFile(reference, "CDS", false);
		final HashMap<String,GffLine> old=unique(before), next=unique(after), gold=unique(truth);
		final HashSet<String> stops=new HashSet<String>(); for(GffLine row : gold.values()){stops.add(stopKey(row));}
		final ByteStreamWriter writer=new ByteStreamWriter(out, false, false, false); writer.start();
		writer.println("change\tcontig\tstart\tstop\tstrand\texact_reference\tstop_reference\tattributes");
		int added=0, removed=0;
		for(GffLine row : before){if(!next.containsKey(key(row))){emit(writer, "removed", row, gold, stops); removed++;}}
		for(GffLine row : after){if(!old.containsKey(key(row))){emit(writer, "added", row, gold, stops); added++;}}
		if(writer.poisonAndWait()){throw new IllegalStateException("Failed writing changed-call evidence");}
		System.out.println("mode\tcalled\treference\texact_tp\texact_precision\texact_recall\tstop_tp\tstop_precision\tstop_recall\tadded\tremoved");
		summary("baseline", old, gold, stops, 0, 0);
		summary("candidate", next, gold, stops, added, removed);
	}

	private static HashMap<String,GffLine> unique(ArrayList<GffLine> rows){
		final HashMap<String,GffLine> result=new HashMap<String,GffLine>();
		for(GffLine row : rows){
			if(row.strand!=0 && row.strand!=1){throw new IllegalArgumentException("CDS comparison requires a declared strand: "+row);}
			if(result.put(key(row), row)!=null){throw new IllegalArgumentException("Duplicate CDS coordinates in comparison: "+key(row));}
		}
		return result;
	}
	private static void summary(String mode, HashMap<String,GffLine> calls, HashMap<String,GffLine> gold,
			HashSet<String> stops, int added, int removed){
		int exact=0, stop=0; final HashSet<String> coveredStops=new HashSet<String>();
		for(GffLine row : calls.values()){
			if(gold.containsKey(key(row))){exact++;}
			final String end=stopKey(row); if(stops.contains(end)){stop++; coveredStops.add(end);}
		}
		assert(exact<=calls.size() && exact<=gold.size() && coveredStops.size()<=stops.size())
			: "Reference agreement counts cannot exceed their distinct-coordinate denominators";
		System.out.println(mode+"\t"+calls.size()+"\t"+gold.size()+"\t"+exact+"\t"+ratio(exact,calls.size())
			+"\t"+ratio(exact,gold.size())+"\t"+stop+"\t"+ratio(stop,calls.size())+"\t"+ratio(coveredStops.size(),stops.size())+"\t"+added+"\t"+removed);
	}
	private static double ratio(int n, int d){return d==0 ? Double.NaN : (double)n/d;}
	private static void emit(ByteStreamWriter writer, String change, GffLine row,
			HashMap<String,GffLine> gold, HashSet<String> stops){
		final ByteBuilder line=new ByteBuilder();
		line.append(change).tab().append(row.seqid).tab().append(row.start).tab().append(row.stop).tab()
			.append(row.strand).tab().append(gold.containsKey(key(row)) ? 1 : 0).tab()
			.append(stops.contains(stopKey(row)) ? 1 : 0).tab().append(row.attributes);
		writer.println(line);
	}
	private static String key(GffLine row){return row.seqid+"\t"+row.strand+"\t"+row.start+"\t"+row.stop;}
	private static String stopKey(GffLine row){return row.seqid+"\t"+row.strand+"\t"+(row.strand==0 ? row.stop : row.start);}
	private ProteinSupportCompare(){}
}
