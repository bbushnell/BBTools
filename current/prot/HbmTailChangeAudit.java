package prot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Joins changed stop-level rows to normalized reference/baseline/candidate rows.
 * Primary stop matches are retained separately from raw exact-coordinate matches.
 * @author Yelan
 */
public final class HbmTailChangeAudit{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace(System.err); System.exit(1);}
	}

	public static void run(String[] args){
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		for(String key : opt.keySet()){
			require(key.equals("stops") || key.equals("reference") || key.equals("baseline")
				|| key.equals("candidate") || key.equals("out"), "Unknown change-audit option: "+key);
		}
		final HashMap<String,Row> ref=loadRows(HmmComparisonData.required(opt, "reference"), "reference");
		final HashMap<String,Row> base=loadRows(HmmComparisonData.required(opt, "baseline"), "baseline");
		final HashMap<String,Row> cand=loadRows(HmmComparisonData.required(opt, "candidate"), "candidate");
		final HashSet<String> refExact=coordinateKeys(ref);
		final ArrayList<Event> events=loadEvents(HmmComparisonData.required(opt, "stops"));
		for(Event e : events){
			validateResolved("reference", resolve(ref, e.referenceIds, "reference", e), e.referencePrimary, e);
			e.baselineRows=resolve(base, e.baselineIds, "baseline", e);
			validateResolved("baseline", e.baselineRows, e.baselinePrimary, e);
			validateResolved("candidate", resolve(cand, e.candidateIds, "candidate", e), e.candidatePrimary, e);
		}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(opt, "out"));
		int rows=0;
		try{
			out.println("event_id\tchange\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\treference_row_ids\tbaseline_row_ids\tcandidate_row_ids\tside\trow_id\texact_raw_reference\tprimary_stop_reference\topposite_lost_exact_overlap_count\topposite_lost_exact_rows\tstart\tend\tlength\tpartial5\tpartial3\tstratum\tnative_score\tcall_bonus\tbbgap_status\tbbgap_rule\tfamily_id\tscore64\tscore64_per_paired\tpaired\tidentity\tterminal20_score64\tbest_other_score64\trunner_supported_score64\tterminal_paired\tquery_suffix\tmodel_suffix\tquery_start\tmodel_start\tlabel_class\tlabel_support\tproduct\tattributes");
			final ByteBuilder bb=new ByteBuilder();
			for(int i=0; i<events.size(); i++){
				final Event e=events.get(i);
				rows+=emit(out, events, i+1, e, "reference", ref, e.referenceIds, refExact, bb);
				rows+=emit(out, events, i+1, e, "baseline", base, e.baselineIds, refExact, bb);
				rows+=emit(out, events, i+1, e, "candidate", cand, e.candidateIds, refExact, bb);
			}
		}finally{HmmComparisonData.close(out);}
		System.err.println("HBM_TAIL_CHANGE_AUDIT_PASS events="+events.size()+" rows="+rows);
	}

	private static int emit(ByteStreamWriter out, ArrayList<Event> events, int eventId, Event e,
			String side, HashMap<String,Row> rows, String ids, HashSet<String> refExact, ByteBuilder bb){
		if(ids.equals(".")){return 0;}
		int emitted=0;
		for(String id : ids.split(",")){
			final Row row=rows.get(id);
			require(row!=null, "Missing "+side+" normalized row "+id+" for "+e.change+" stop "+e.stopKey());
			require(row.matches(e), side+" normalized row "+id+" does not match event "+e.stopKey());
			final String overlap=side.equals("candidate") ? oppositeLost(events, row, refExact) : ".";
			bb.clear().append(eventId).tab().append(e.change).tab().append(e.contig).tab().append(e.strand).tab().append(e.terminal).tab()
				.append(e.referenceRaw).tab().append(e.referencePrimary).tab().append(e.baselineRaw).tab().append(e.baselinePrimary).tab()
				.append(e.candidateRaw).tab().append(e.candidatePrimary).tab()
				.append(e.referenceIds).tab().append(e.baselineIds).tab().append(e.candidateIds).tab()
				.append(side).tab().append(row.id).tab().append(exact(row, refExact) ? 1 : 0).tab().append(e.referencePrimary>0 ? 1 : 0).tab()
				.append(overlap.equals(".") ? 0 : overlap.split(",").length).tab().append(overlap).tab()
				.append(row.start).tab().append(row.end).tab().append(row.end-row.start+1).tab()
				.append(row.partial5).tab().append(row.partial3).tab().append(row.stratum).tab()
				.append(attr(row, "bbgap_native_score")).tab().append(attr(row, "bbgap_call_bonus")).tab()
				.append(attr(row, "bbgap_status")).tab().append(attr(row, "bbgap_rule")).tab()
				.append(attr(row, "bbgap_family_id")).tab().append(attr(row, "bbgap_score64")).tab()
				.append(scorePerPaired(row)).tab().append(attr(row, "bbgap_tail_paired")).tab().append(attr(row, "bbgap_identity")).tab()
				.append(attr(row, "bbgap_tail_terminal_score64")).tab().append(attr(row, "bbgap_best_other_score64")).tab()
				.append(attr(row, "bbgap_runner_supported_score64")).tab().append(attr(row, "bbgap_tail_terminal_paired")).tab()
				.append(attr(row, "bbgap_tail_query_suffix")).tab().append(attr(row, "bbgap_tail_model_suffix")).tab()
				.append(attr(row, "bbgap_tail_query_start")).tab().append(attr(row, "bbgap_tail_model_start")).tab()
				.append(attr(row, "bbgap_label_class")).tab().append(attr(row, "bbgap_label_support")).tab()
				.append(attr(row, "bbgap_product")).tab().append(row.attributes);
			out.println(bb);
			emitted++;
		}
		return emitted;
	}

	private static HashMap<String,Row> loadRows(String path, String side){
		final HashMap<String,Row> map=new HashMap<String,Row>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			byte[] line=bf.nextLine();
			require(line!=null, "Empty normalized row table: "+path);
			lp.set(line);
			require(equals(line, ROW_HEADER), "Expected normalized stops-v1 row schema in "+path);
			for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==13, "Malformed "+side+" row in "+path+": "+new String(line));
				require(lp.termEquals("stops-v1", 0), "Malformed "+side+" schema in "+path+": "+new String(line));
				final Row row=new Row(lp.parseString(1), lp.parseString(3), lp.parseInt(4),
					lp.parseInt(5), lp.parseInt(6), lp.parseInt(9), lp.parseInt(10),
					lp.parseString(11), lp.parseString(12));
				String key=row.id;
				if(map.containsKey(key)){
					require(key.equals("0"), "Duplicate "+side+" normalized input_row: "+key);
					key="0:"+map.size();
				}
				if(key.equals("0")){key="0:"+map.size();}
				map.put(key, row);
			}
		}finally{require(!bf.close(), "Failed reading normalized row table: "+path);}
		return map;
	}

	private static ArrayList<Event> loadEvents(String path){
		final ArrayList<Event> list=new ArrayList<Event>();
		final HashSet<String> keys=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			byte[] line=bf.nextLine();
			require(line!=null, "Empty stop table: "+path);
			lp.set(line);
			require(equals(line, STOP_HEADER), "Expected paired stops-v1 event schema in "+path);
			for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==14, "Malformed stop row in "+path+": "+new String(line));
				require(lp.termEquals("stops-v1", 0), "Malformed stop schema in "+path+": "+new String(line));
				final String key=lp.parseString(1)+"\t"+lp.parseInt(2)+"\t"+lp.parseInt(3);
				require(keys.add(key), "Duplicate stop event key: "+key);
				final String change=lp.parseString(10);
				if(change.equals("unchanged") || change.equals("excluded_or_uncalled")){continue;}
				require(change.equals("gained") || change.equals("lost") || change.equals("same_stop_changed_starts"),
					"Unknown paired stop change: "+change);
				final Event e=new Event(lp.parseString(1), lp.parseInt(2), lp.parseInt(3),
					lp.parseInt(4), lp.parseInt(5), lp.parseInt(6), lp.parseInt(7),
					lp.parseInt(8), lp.parseInt(9), change,
					lp.parseString(11), lp.parseString(12), lp.parseString(13));
				validateIds("reference", e.referenceIds, e.referenceRaw, e);
				validateIds("baseline", e.baselineIds, e.baselineRaw, e);
				validateIds("candidate", e.candidateIds, e.candidateRaw, e);
				list.add(e);
			}
		}finally{require(!bf.close(), "Failed reading stop table: "+path);}
		return list;
	}

	private static HashSet<String> coordinateKeys(HashMap<String,Row> rows){
		final HashSet<String> set=new HashSet<String>();
		for(Row row : rows.values()){set.add(row.coordinateKey());}
		return set;
	}
	private static boolean exact(Row row, HashSet<String> refExact){
		return refExact.contains(row.coordinateKey());
	}
	private static String oppositeLost(ArrayList<Event> events, Row row, HashSet<String> refExact){
		final ByteBuilder bb=new ByteBuilder();
		for(Event e : events){
			if(!e.change.equals("lost") || e.baselineRows==null){continue;}
			for(Row old : e.baselineRows){
				if(old.primary() && row.strand!=old.strand && row.contig.equals(old.contig) && row.start<=old.end && row.end>=old.start && exact(old, refExact)){
					if(bb.length()>0){bb.append(',');}
					bb.append(old.id);
				}
			}
		}
		return bb.length()<1 ? "." : bb.toString();
	}
	private static ArrayList<Row> resolve(HashMap<String,Row> rows, String ids, String side, Event e){
		final ArrayList<Row> list=new ArrayList<Row>();
		if(ids.equals(".")){return list;}
		for(String id : ids.split(",")){
			final Row row=rows.get(id);
			require(row!=null, "Missing "+side+" normalized row "+id+" for "+e.change+" stop "+e.stopKey());
			require(row.matches(e), side+" normalized row "+id+" does not match event "+e.stopKey());
			list.add(row);
		}
		return list;
	}
	private static void validateIds(String side, String ids, int raw, Event e){
		final int count=countIds(ids);
		require(count==raw, side+" row-count mismatch for "+e.change+" stop "+e.stopKey());
		if(ids.equals(".")){return;}
		final HashSet<String> seen=new HashSet<String>();
		for(String id : ids.split(",")){
			require(!id.equals("0") && id.length()>0, "Ambiguous "+side+" normalized row ID for "+e.change+" stop "+e.stopKey());
			require(seen.add(id), "Duplicate "+side+" normalized row ID "+id+" for "+e.change+" stop "+e.stopKey());
		}
	}
	private static void validateResolved(String side, ArrayList<Row> rows, int primary, Event e){
		int count=0;
		for(Row row : rows){if(row.primary()){count++;}}
		require(count==primary, side+" primary-count mismatch for "+e.change+" stop "+e.stopKey());
	}
	private static int countIds(String ids){return ids.equals(".") ? 0 : ids.split(",").length;}
	private static boolean equals(byte[] line, String expected){
		return line.length==expected.length() && expected.contentEquals(new String(line));
	}
	private static String attr(Row row, String key){
		final String a=row.attributes;
		for(int start=0; start<a.length(); ){
			int end=a.indexOf(';', start); if(end<0){end=a.length();}
			final int eq=a.indexOf('=', start);
			if(eq>start && eq<end && a.regionMatches(start, key, 0, key.length()) && eq==start+key.length()){
				return a.substring(eq+1, end);
			}
			start=end+1;
		}
		return ".";
	}
	private static String scorePerPaired(Row row){
		final String score=attr(row, "bbgap_score64"), paired=attr(row, "bbgap_tail_paired");
		if(score.equals(".") || paired.equals(".")){return ".";}
		final int p=Integer.parseInt(paired);
		return p<1 ? "." : Float.toString(Integer.parseInt(score)/(float)p);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	private static final class Event{
		Event(String contig_, int strand_, int terminal_, int rr_, int rp_, int br_, int bp_,
				int cr_, int cp_, String change_, String ref_, String base_, String cand_){
			contig=contig_; strand=strand_; terminal=terminal_; referenceRaw=rr_; referencePrimary=rp_;
			baselineRaw=br_; baselinePrimary=bp_; candidateRaw=cr_; candidatePrimary=cp_;
			change=change_; referenceIds=ref_; baselineIds=base_; candidateIds=cand_;
		}
		String stopKey(){return contig+"\t"+strand+"\t"+terminal;}
		final String contig, change, referenceIds, baselineIds, candidateIds;
		final int strand, terminal, referenceRaw, referencePrimary, baselineRaw, baselinePrimary, candidateRaw, candidatePrimary;
		ArrayList<Row> baselineRows;
	}
	private static final class Row{
		Row(String id_, String contig_, int strand_, int start_, int end_, int partial5_, int partial3_,
				String stratum_, String attributes_){
			id=id_; contig=contig_; strand=strand_; start=start_; end=end_;
			partial5=partial5_; partial3=partial3_; stratum=stratum_; attributes=attributes_;
		}
		boolean matches(Event e){return strand==e.strand && terminal()==e.terminal && contig.equals(e.contig);}
		boolean primary(){return stratum.equals("primary");}
		int terminal(){return strand==0 ? end : start;}
		String coordinateKey(){return contig+"\t"+strand+"\t"+start+"\t"+end;}
		final String id, contig, stratum, attributes;
		final int strand, start, end, partial5, partial3;
	}
	private static final String ROW_HEADER="schema\tinput_row\ttype\tcontig\tstrand\tstart\tend\tphase\ttable\tpartial5\tpartial3\tstratum\tattributes";
	private static final String STOP_HEADER="schema\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tchange\treference_rows\tbaseline_rows\tcandidate_rows";
	private HbmTailChangeAudit(){}
}
