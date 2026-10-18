package prot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.GffLine;
import parse.LineParser1;
import structures.ByteBuilder;

/** Joins failed tail counterfactual rows to raw reference and selected-call GFF stops.
 * The reference columns are raw coordinate concordance, not normalized truth.
 * @author Yelan
 */
public final class HbmTailLossAudit{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace(System.err); System.exit(1);}
	}

	private static void run(String[] args){
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		for(String key : opt.keySet()){require(key.equals("list") || key.equals("out"), "Unknown loss-audit option: "+key);}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(opt, "out"));
		out.println("genome\tquery\tfamily_id\tcontig\tstart\tend\tstrand\tin_final_call\tsame_final_stop\traw_ref_exact\traw_ref_stop\traw_ref_stop_count\tscore64\tidentity\tpaired_fraction\tcore_coverage\tpaired\tpaired_core\tterminal_paired\tterminal_identical\tquery_terminal_paired\tquery_terminal_identical\tquery_suffix\treference_suffix\tquery_start\treference_start\tshortlist_rank\tterminal20_score64\tterminal40_score64\tpaired_identity");
		long files=0, losses=0;
		for(String[] row : HmmComparisonData.rows(HmmComparisonData.required(opt, "list"))){
			if(row.length>0 && row[0].equals("genome")){continue;}
			require(row.length==5, "Loss-audit list rows must have genome, rows, support, calls, ref");
			losses+=audit(row[0], row[1], row[2], row[3], row[4], out);
			files++;
		}
		HmmComparisonData.close(out);
		require(files>0 && losses>0, "Loss audit requires at least one input and one loss");
		System.err.println("HBM_TAIL_LOSS_AUDIT_PASS genomes="+files+" losses="+losses+" raw_reference=true");
	}

	private static int audit(String genome, String rowsPath, String supportPath, String callsPath, String refPath, ByteStreamWriter out){
		final HashMap<String,Support> support=loadSupport(genome, supportPath);
		final GffIndex calls=new GffIndex(GffLine.loadGffFile(callsPath, "CDS", false));
		final GffIndex refs=new GffIndex(GffLine.loadGffFile(refPath, "CDS", false));
		final ByteFile bf=ByteFile.makeByteFile(rowsPath, false);
		final LineParser1 lp=new LineParser1('\t');
		int losses=0;
		try{
			byte[] line=bf.nextLine();
			require(line!=null, "Empty tail diagnostic table: "+rowsPath);
			lp.set(line);
			require(lp.terms()==27 && lp.termEquals("query", 0) && lp.termEquals("family_id", 2)
				&& lp.termEquals("variant", 3) && lp.termEquals("status", 5)
				&& lp.termEquals("terminal_paired", 12) && lp.termEquals("terminal_identical", 13)
				&& lp.termEquals("shortlist_rank", 22) && lp.termEquals("terminal20_score64", 23)
				&& lp.termEquals("terminal40_score64", 24) && lp.termEquals("paired_identity", 25)
				&& lp.termEquals("counterfactual_pass", 26), "Expected HbmTailRobustness 27-column schema");
			for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==27, "Malformed tail row in "+rowsPath+": "+new String(line));
				if(!lp.termEquals("original", 3) || lp.termEquals("true", 26)){continue;}
				require(lp.termEquals("false", 26), "Counterfactual flag must be true or false");
				final String query=lp.parseString(0);
				require(lp.termEquals("stop_supported", 5), "Loss audit expects original current-supported rows: "+query);
				final Support s=support.remove(query);
				require(s!=null, "Tail loss missing from support map: "+query);
				require(s.family==lp.parseInt(2), "Family mismatch for "+query);
				final ByteBuilder bb=new ByteBuilder();
				bb.append(genome).tab().append(query).tab().append(s.family).tab().append(s.contig).tab()
					.append(s.start).tab().append(s.stop).tab().append(s.strand).tab()
					.append(calls.exact.contains(s.exactKey()) ? 1 : 0).tab()
					.append(calls.stops.containsKey(s.stopKey()) ? 1 : 0).tab()
					.append(refs.exact.contains(s.exactKey()) ? 1 : 0).tab()
					.append(refs.stops.containsKey(s.stopKey()) ? 1 : 0).tab()
					.append(refs.stops.getOrDefault(s.stopKey(), 0));
				for(int i=6; i<=19; i++){bb.tab().append(lp.parseString(i));}
				bb.tab().append(lp.parseString(22)).tab().append(lp.parseString(23))
					.tab().append(lp.parseString(24)).tab().append(lp.parseString(25));
				out.println(bb);
				losses++;
			}
		}finally{require(!bf.close(), "Failed reading tail diagnostic table: "+rowsPath);}
		return losses;
	}

	private static HashMap<String,Support> loadSupport(String genome, String path){
		final HashMap<String,Support> map=new HashMap<String,Support>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			byte[] line=bf.nextLine();
			require(line!=null, "Empty current-support table: "+path);
			lp.set(line);
			require(lp.terms()==18 && lp.termEquals("query", 0) && lp.termEquals("genome", 1)
				&& lp.termEquals("contig", 2) && lp.termEquals("start", 3) && lp.termEquals("stop", 4)
				&& lp.termEquals("strand", 5) && lp.termEquals("family_id", 6),
				"Expected CallGenes current-support 18-column schema");
			for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==18, "Malformed current-support row in "+path+": "+new String(line));
				require(lp.termEquals(genome, 1), "Genome mismatch in support table: "+lp.parseString(1));
				final Support s=new Support(lp.parseString(2), lp.parseInt(3), lp.parseInt(4),
					strand(lp.parseString(5)), lp.parseInt(6));
				require(map.put(lp.parseString(0), s)==null, "Duplicate support query: "+lp.parseString(0));
			}
		}finally{require(!bf.close(), "Failed reading current-support table: "+path);}
		return map;
	}

	private static int strand(String s){
		if(s.equals("+")){return GffLine.PLUS;}
		if(s.equals("-")){return GffLine.MINUS;}
		throw new IllegalArgumentException("Expected +/- strand: "+s);
	}

	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	private static final class GffIndex{
		GffIndex(ArrayList<GffLine> rows){
			for(GffLine row : rows){
				if(row.strand!=GffLine.PLUS && row.strand!=GffLine.MINUS){
					throw new IllegalArgumentException("CDS loss audit requires a declared strand: "+row);
				}
				exact.add(key(row));
				final String stop=stop(row);
				stops.put(stop, stops.getOrDefault(stop, 0)+1);
			}
		}
		final HashSet<String> exact=new HashSet<String>();
		final HashMap<String,Integer> stops=new HashMap<String,Integer>();
	}
	private static final class Support{
		Support(String contig_, int start_, int stop_, int strand_, int family_){
			contig=contig_; start=start_; stop=stop_; strand=strand_; family=family_;
		}
		String exactKey(){return contig+"\t"+strand+"\t"+start+"\t"+stop;}
		String stopKey(){return contig+"\t"+strand+"\t"+(strand==GffLine.PLUS ? stop : start);}
		final String contig;
		final int start, stop, strand, family;
	}
	private static String key(GffLine row){return row.seqid+"\t"+row.strand+"\t"+row.start+"\t"+row.stop;}
	private static String stop(GffLine row){return row.seqid+"\t"+row.strand+"\t"+(row.strand==GffLine.PLUS ? row.stop : row.start);}
	private HbmTailLossAudit(){}
}
