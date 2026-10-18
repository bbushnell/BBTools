package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import shared.Tools;
import structures.ByteBuilder;

/** Joins supported-stop fate diagnostics to normalized selected/reference rows.
 * Observational only: a selected overlap is not a DP causality claim.
 * mode=changes joins evaluator stop gains/losses, including keys absent from supported-only traces.
 * @author Yelan, Sayu
 */
public final class ProteinSupportFateReport{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace(System.err); System.exit(1);}
	}

	public static void run(String[] args){
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		for(String key : opt.keySet()){
			require(key.equals("manifest") || key.equals("out") || key.equals("mode"), "Unknown fate-report option: "+key);
		}
		final String mode=opt.containsKey("mode") ? opt.get("mode") : "fates";
		require(mode.equals("fates") || mode.equals("changes"), "Fate-report mode must be fates or changes");
		if(mode.equals("changes")){writeChanges(opt); return;}
		final ArrayList<Run> runs=loadManifest(HmmComparisonData.required(opt, "manifest"));
		final String out=HmmComparisonData.required(opt, "out");
		final String[] outputs={out+".stops.tsv", out+".losses.tsv", out+".summary.tsv"};
		require(Tools.testOutputFiles(false, false, false, outputs), "Fate-report outputs must be new writable files");
		write(runs, outputs);
		int fateRows=0, omitted=0, losses=0;
		for(Run run : runs){
			fateRows+=run.fates.size(); omitted+=run.omitted; losses+=run.referenceLosses;
		}
		System.err.println("PROTEIN_SUPPORT_FATE_REPORT_PASS genomes="+runs.size()+" supported_stops="+fateRows
			+" omitted_supported_stops="+omitted+" reference_losses="+losses+" selected_keys_conserved=true");
	}

	private static ArrayList<Run> loadManifest(String path){
		final ArrayList<Run> runs=new ArrayList<Run>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(bf, lp, "genome", "fate", "reference", "candidate");
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==4, "Malformed fate manifest row in "+path+": "+new String(line));
				final Run run=new Run(lp.parseString(0));
				require(!run.genome.isEmpty() && seen.add(run.genome), "Duplicate/empty genome in fate manifest: "+run.genome);
				run.fates=loadFates(lp.parseString(1));
				run.references=loadReference(lp.parseString(2));
				run.selected=loadSelected(lp.parseString(3));
				validate(run);
				runs.add(run);
			}
		}finally{require(!bf.close(), "Failed reading fate manifest: "+path);}
		require(!runs.isEmpty(), "Empty fate manifest: "+path);
		return runs;
	}

	private static ArrayList<Fate> loadFates(String path){
		final ArrayList<Fate> rows=new ArrayList<Fate>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(bf, lp, FATE_HEADER);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==FATE_HEADER.length && lp.termEquals("support-fate-v1", 0),
					"Malformed support-fate row in "+path+": "+new String(line));
				final Fate fate=new Fate();
				final String oldKey=lp.parseString(1);
				fate.contig=unesc(lp.parseString(2)); fate.strand=lp.parseInt(3); fate.terminal=lp.parseInt(4);
				fate.key=esc(fate.contig)+":"+fate.strand+":"+fate.terminal;
				fate.longestStart=lp.parseString(5); fate.longestEnd=lp.parseString(6);
				fate.family=lp.parseString(9); fate.identity=lp.parseString(11); fate.score64=lp.parseString(12);
				fate.admitted=lp.parseInt(17); fate.rejectedStart=lp.parseInt(18); fate.rejectedLength=lp.parseInt(19);
				fate.rejectedLimit=lp.parseInt(20); fate.softRejected=lp.parseInt(21); fate.retained=lp.parseInt(22);
				fate.selected=lp.parseInt(23); fate.fate=lp.parseString(24);
				fate.predicates=lp.parseString(25); fate.bestStart=lp.parseString(26); fate.bestEnd=lp.parseString(27);
				fate.bestRawStart=lp.parseString(28); fate.bestRawStop=lp.parseString(29); fate.bestInner=lp.parseString(30);
				fate.bestNative=lp.parseString(31); fate.bestBonus=lp.parseString(32); fate.bestTotal=lp.parseString(33);
				fate.overlapBasis=lp.parseString(34); fate.overlapStart=lp.parseString(35); fate.overlapEnd=lp.parseString(36);
				fate.overlaps=lp.parseString(39);
				require(fate.strand==0 || fate.strand==1, "Invalid fate strand: "+fate.key);
				require(fate.key.equals(oldKey), "stop_id must match decoded contig, strand and terminal: "+oldKey);
				require(fate.terminal>0 && fate.admitted>=0 && fate.rejectedStart>=0 && fate.rejectedLength>=0
					&& fate.rejectedLimit>=0 && fate.softRejected>=0 && fate.retained>=0
					&& (fate.selected==0 || fate.selected==1) && fate.selected<=fate.retained, "Invalid fate counts: "+fate.key);
				require(fate.admitted==fate.softRejected+fate.retained, "Fate rejected+retained must equal admitted: "+fate.key);
				require(predicateSum(fate.predicates)==fate.softRejected,
					"Pruning predicates must sum to soft-pruned candidates: "+fate.key);
				require(seen.add(fate.key), "Duplicate fate stop: "+fate.key);
				rows.add(fate);
			}
		}finally{require(!bf.close(), "Failed reading fate rows: "+path);}
		Collections.sort(rows, ORDER);
		return rows;
	}

	private static HashMap<String,Reference> loadReference(String path){
		final HashMap<String,Reference> refs=new HashMap<String,Reference>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(bf, lp, ROW_HEADER);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==ROW_HEADER.length && lp.termEquals("stops-v1", 0),
					"Malformed reference row in "+path+": "+new String(line));
				if(!lp.termEquals("CDS", 2)){continue;}
				final int strand=lp.parseInt(4);
				if(strand!=0 && strand!=1){continue;}
				final String key=key(lp);
				Reference ref=refs.get(key);
				if(ref==null){ref=new Reference(); refs.put(key, ref);}
				ref.raw++;
				if(lp.termEquals("primary", 11)){ref.primary++;}
				final String stratum=lp.parseString(11);
				final Integer old=ref.strata.get(stratum);
				ref.strata.put(stratum, old==null ? 1 : old+1);
			}
		}finally{require(!bf.close(), "Failed reading reference rows: "+path);}
		return refs;
	}

	private static HashSet<String> loadSelected(String path){
		final HashSet<String> selected=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(bf, lp, ROW_HEADER);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==ROW_HEADER.length && lp.termEquals("stops-v1", 0),
					"Malformed candidate row in "+path+": "+new String(line));
				if(lp.termEquals("CDS", 2) && attr(lp.parseString(12), "bbgap_status").equals("stop_supported")){
					final String key=key(lp);
					require(selected.add(key), "Duplicate selected stop-supported key: "+key);
				}
			}
		}finally{require(!bf.close(), "Failed reading candidate rows: "+path);}
		return selected;
	}

	private static void validate(Run run){
		final HashSet<String> selected=new HashSet<String>(run.selected);
		for(Fate fate : run.fates){
			final boolean chosen=selected.remove(fate.key);
			require((fate.selected>0)==chosen, "Fate/GFF selected-key mismatch for "+run.genome+"/"+fate.key);
			final boolean named=fate.fate.equals("selected") || fate.fate.equals("dp_not_selected")
				|| fate.fate.equals("no_starts") || fate.fate.equals("soft_pruned");
			require(named, "Unknown fate class: "+fate.fate);
			final String expected=(fate.selected>0 ? "selected" : fate.retained>0 ? "dp_not_selected"
				: fate.admitted==0 ? "no_starts" : "soft_pruned");
			require(fate.fate.equals(expected), "Fate class disagrees with counts for "+run.genome+"/"+fate.key);
			if(fate.selected==0){
				run.omitted++;
				if(reference(run, fate).primary>0){run.referenceLosses++;}
			}
		}
		require(selected.isEmpty(), "Selected stop-supported rows missing from fate output for "+run.genome+": "+selected);
	}

	private static void write(ArrayList<Run> runs, String[] outputs){
		final ByteStreamWriter stops=HmmComparisonData.writer(outputs[0]);
		final ByteStreamWriter losses=HmmComparisonData.writer(outputs[1]);
		final ByteStreamWriter summary=HmmComparisonData.writer(outputs[2]);
		final ByteBuilder row=new ByteBuilder();
		try{
			stops.println(STOP_HEADER);
			losses.println(STOP_HEADER);
			summary.println("genome\tsupported_stops\tselected_supported_stops\tomitted_supported_stops\treference_primary_supported\treference_primary_selected\treference_primary_not_selected\tno_start_stops\tsoft_pruned_stops\tdp_not_selected_stops");
			for(Run run : runs){
				int selected=0, refPrimary=0, refSelected=0, noStarts=0, softPruned=0, dp=0;
				for(Fate fate : run.fates){
					final Reference ref=reference(run, fate);
					if(fate.selected>0){selected++;}
					if(ref.primary>0){refPrimary++; if(fate.selected>0){refSelected++;}}
					if(fate.fate.equals("no_starts")){noStarts++;}
					if(fate.fate.equals("soft_pruned")){softPruned++;}
					if(fate.fate.equals("dp_not_selected")){dp++;}
					row.clear().append(run.genome).tab().append(fate.key).tab().append(esc(fate.contig)).tab().append(fate.strand)
						.tab().append(fate.terminal).tab().append(ref.raw).tab().append(ref.primary).tab().append(ref.state())
						.tab().append(ref.strata()).tab().append(fate.admitted).tab().append(fate.rejectedStart)
						.tab().append(fate.rejectedLength).tab().append(fate.rejectedLimit).tab().append(fate.softRejected)
						.tab().append(fate.retained).tab().append(fate.selected).tab().append(fate.fate)
						.tab().append(fate.predicates).tab().append(fate.family).tab().append(fate.identity)
						.tab().append(fate.score64).tab().append(fate.longestStart).tab().append(fate.longestEnd)
						.tab().append(fate.bestStart).tab().append(fate.bestEnd).tab().append(fate.bestRawStart)
						.tab().append(fate.bestRawStop).tab().append(fate.bestInner).tab().append(fate.bestNative)
						.tab().append(fate.bestBonus).tab().append(fate.bestTotal).tab().append(fate.overlapBasis)
						.tab().append(fate.overlapStart).tab().append(fate.overlapEnd).tab().append(fate.overlaps);
					stops.println(row);
					if(fate.selected==0){losses.println(row);}
				}
				summary.println(row.clear().append(run.genome).tab().append(run.fates.size()).tab().append(selected)
					.tab().append(run.omitted).tab().append(refPrimary).tab().append(refSelected)
					.tab().append(refPrimary-refSelected).tab().append(noStarts).tab().append(softPruned).tab().append(dp));
			}
		}finally{
			HmmComparisonData.close(stops);
			HmmComparisonData.close(losses);
			HmmComparisonData.close(summary);
		}
	}

	private static Reference reference(Run run, Fate fate){
		final Reference ref=run.references.get(fate.key);
		return ref==null ? Reference.EMPTY : ref;
	}

	/** Separate opt-in outputs; legacy manifest, joins and output bytes are unchanged. */
	private static void writeChanges(HashMap<String,String> opt){
		final String path=HmmComparisonData.required(opt, "manifest"), out=HmmComparisonData.required(opt, "out");
		final ArrayList<ChangeRun> runs=new ArrayList<ChangeRun>();
		final HashSet<String> genomes=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(bf, lp, "genome", "fate", "comparison");
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==3, "Malformed changes manifest row in "+path);
				final String genome=lp.parseString(0);
				require(!genome.isEmpty() && genomes.add(genome), "Duplicate/empty genome in changes manifest: "+genome);
				runs.add(loadChanges(genome, lp.parseString(1), lp.parseString(2)));
			}
		}finally{require(!bf.close(), "Failed reading changes manifest: "+path);}
		require(!runs.isEmpty(), "Empty changes manifest: "+path);
		require(Tools.testOutputFiles(false, false, false, out+".changes.tsv", out+".summary.tsv"),
			"Change-report outputs must be new writable files");
		final ByteStreamWriter changes=HmmComparisonData.writer(out+".changes.tsv");
		final ByteStreamWriter summary=HmmComparisonData.writer(out+".summary.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			changes.println("genome\tstop_id\tcontig\tstrand\tterminal\tchange\treference_raw_rows\treference_primary_rows"
				+"\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tfate_traced"
				+"\tselected_supported\tfate\tfirst_pruning_predicates\tbest_inner\tbest_native\tbest_bonus");
			summary.println("genome\tchange\treference_stops\tfate_supported_stops\tselected_supported_stops"
				+"\tno_starts\tsoft_pruned\tdp_not_selected\tuntraced");
			for(ChangeRun run : runs){
				final long[][] counts=new long[2][5];
				for(Change change : run.changes.values()){
					final Fate fate=change.fate;
					final int category=fateCategory(fate);
					counts[change.lost ? 0 : 1][category]++;
					row.clear().append(run.genome).tab().append(change.key).tab().append(esc(change.contig))
						.tab().append(change.strand).tab().append(change.terminal).tab().append(change.lost ? "lost" : "gained");
					for(int count : change.rows){row.tab().append(count);}
					row.tab().append(fate==null ? 0 : 1).tab().append(fate==null ? "." : Integer.toString(fate.selected))
						.tab().append(fate==null ? "untraced" : fate.fate).tab().append(fate==null ? "." : fate.predicates)
						.tab().append(fate==null ? "." : fate.bestInner).tab().append(fate==null ? "." : fate.bestNative)
						.tab().append(fate==null ? "." : fate.bestBonus);
					changes.println(row);
				}
				long all=0;
				for(int arm=0; arm<counts.length; arm++){
					final long[] c=counts[arm];
					final long traced=c[0]+c[1]+c[2]+c[3], total=traced+c[4];
					all+=total;
					summary.println(row.clear().append(run.genome).tab().append(arm==0 ? "lost" : "gained")
						.tab().append(total).tab().append(traced).tab().append(c[0]).tab().append(c[1])
						.tab().append(c[2]).tab().append(c[3]).tab().append(c[4]));
				}
				assert(all==run.changes.size()) : "Every reference change must belong to exactly one fate partition: "+run.genome;
			}
		}finally{HmmComparisonData.close(changes); HmmComparisonData.close(summary);}
		System.err.println("PROTEIN_SUPPORT_FATE_CHANGES_PASS genomes="+runs.size()+" scope=primary_reference_stop_keys untraced_cause=unknown");
	}

	private static ChangeRun loadChanges(String genome, String fatePath, String comparisonPath){
		final ChangeRun run=new ChangeRun(genome);
		final HashMap<String,Fate> fates=new HashMap<String,Fate>();
		for(Fate fate : loadFates(fatePath)){fateCategory(fate); fates.put(fate.key, fate);}
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(comparisonPath, false);
		final LineParser1 lp=new LineParser1('\t');
		final int[] rows=new int[6];
		try{
			header(bf, lp, COMPARISON_HEADER);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line);
				require(lp.terms()==COMPARISON_HEADER.length && lp.termEquals("stops-v1", 0), "Malformed stop comparison row: "+comparisonPath);
				final String contig=lp.parseString(1);
				final int strand=lp.parseInt(2), terminal=lp.parseInt(3);
				require(!contig.isEmpty() && terminal>0, "Invalid comparison stop coordinates");
				final String key=esc(contig)+":"+strand+":"+terminal;
				require(seen.add(key), "Duplicate comparison stop: "+key);
				for(int i=0; i<rows.length; i++){
					rows[i]=lp.parseInt(4+i);
					require(rows[i]>=0 && (i%2==0 || rows[i]<=rows[i-1]), "Invalid raw/primary comparison counts: "+key);
				}
				final String change=lp.parseString(10);
				final boolean before=rows[3]>0, after=rows[5]>0;
				final boolean valid=before ? (after ? change.equals("unchanged") || change.equals("same_stop_changed_starts")
					: change.equals("lost")) : after ? change.equals("gained") : change.equals("excluded_or_uncalled");
				require(valid, "Change class disagrees with primary counts: "+key);
				if(strand!=0 && strand!=1){
					// ProteinStopEvaluation.classify excludes unknown strands, but Dataset.index
					// and writeStops retain their raw endpoints. These cannot have supported fates.
					require(rows[1]==0 && rows[3]==0 && rows[5]==0, "Primary comparison rows require a known strand: "+key);
					continue;
				}
				final Fate fate=fates.get(key);
				if(fate!=null){require((fate.selected>0)==(rows[4]>0), "Fate/comparison candidate-presence mismatch: "+key);}
				if(rows[1]>0 && before!=after){
					final Change c=new Change();
					c.key=key; c.contig=contig; c.strand=strand; c.terminal=terminal;
					c.lost=before; c.rows=rows.clone(); c.fate=fate;
					run.changes.put(key, c);
				}
			}
		}finally{require(!bf.close(), "Failed reading comparison rows: "+comparisonPath);}
		for(Fate fate : fates.values()){
			require(fate.selected==0 || seen.contains(fate.key), "Selected fate stop missing from comparison: "+fate.key);
		}
		return run;
	}

	/** Category is observed candidate fate, never a claim about why an untraced key disappeared. */
	private static int fateCategory(Fate fate){
		if(fate==null){return 4;}
		final int category=fate.selected>0 ? 0 : fate.retained>0 ? 3 : fate.admitted==0 ? 1 : 2;
		final String expected=category==0 ? "selected" : category==1 ? "no_starts" : category==2 ? "soft_pruned" : "dp_not_selected";
		require(fate.fate.equals(expected), "Fate class disagrees with counts: "+fate.key);
		return category;
	}

	private static String key(LineParser1 lp){
		final int strand=lp.parseInt(4), terminal=(strand==0 ? lp.parseInt(6) : lp.parseInt(5));
		require(strand==0 || strand==1, "Only known strands have stop keys");
		return esc(lp.parseString(3))+":"+strand+":"+terminal;
	}

	private static String attr(String a, String key){
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

	private static int predicateSum(String predicates){
		if(predicates.equals(".")){return 0;}
		int sum=0;
		final HashSet<String> seen=new HashSet<String>();
		for(int start=0; start<predicates.length(); ){
			int end=predicates.indexOf(',', start);
			if(end<0){end=predicates.length();}
			final int colon=predicates.lastIndexOf(':', end-1);
			require(colon>start && colon<end-1, "Malformed pruning predicate token: "+predicates);
			require(seen.add(predicates.substring(start, colon)), "Duplicate pruning predicate token: "+predicates);
			final int count=Integer.parseInt(predicates.substring(colon+1, end));
			require(count>0, "Nonpositive pruning predicate count: "+predicates);
			sum+=count;
			start=end+1;
		}
		return sum;
	}

	private static String esc(String s){
		return s.replace("%", "%25").replace(":", "%3A").replace(",", "%2C")
			.replace("\t", "%09").replace("\n", "%0A").replace("\r", "%0D");
	}
	private static String unesc(String s){
		final StringBuilder sb=new StringBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char c=s.charAt(i);
			if(c!='%'){sb.append(c); continue;}
			require(i+2<s.length(), "Incomplete percent escape in diagnostic contig: "+s);
			final String code=s.substring(i+1, i+3).toUpperCase(java.util.Locale.ROOT);
			if(code.equals("25")){sb.append('%');}
			else if(code.equals("3A")){sb.append(':');}
			else if(code.equals("2C")){sb.append(',');}
			else if(code.equals("09")){sb.append('\t');}
			else if(code.equals("0A")){sb.append('\n');}
			else if(code.equals("0D")){sb.append('\r');}
			else{throw new IllegalArgumentException("Unknown percent escape in diagnostic contig: %"+code);}
			i+=2;
		}
		return sb.toString();
	}

	private static void header(ByteFile bf, LineParser1 lp, String... names){
		final byte[] line=bf.nextLine();
		require(line!=null, "Missing TSV header");
		lp.set(line);
		require(lp.terms()==names.length, "Unexpected TSV width");
		for(int i=0; i<names.length; i++){require(lp.termEquals(names[i], i), "Unexpected header column "+i);}
	}

	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	private static final class Run{
		Run(String genome_){genome=genome_;}
		final String genome;
		ArrayList<Fate> fates;
		HashMap<String,Reference> references;
		HashSet<String> selected;
		int omitted, referenceLosses;
	}

	private static final class Fate{
		String key, contig, longestStart, longestEnd, family, identity, score64;
		String fate, predicates, bestStart, bestEnd, bestRawStart, bestRawStop, bestInner, bestNative, bestBonus, bestTotal;
		String overlapBasis, overlapStart, overlapEnd, overlaps;
		int strand, terminal, admitted, rejectedStart, rejectedLength, rejectedLimit, softRejected, retained, selected;
	}

	private static final class ChangeRun{
		ChangeRun(String genome_){genome=genome_;}
		final String genome;
		final TreeMap<String,Change> changes=new TreeMap<String,Change>();
	}
	private static final class Change{
		String key, contig;
		int strand, terminal;
		boolean lost;
		int[] rows;
		Fate fate;
	}

	private static final class Reference{
		static final Reference EMPTY=new Reference();
		int raw, primary;
		final TreeMap<String,Integer> strata=new TreeMap<String,Integer>();
		String state(){return primary>0 ? "primary" : raw>0 ? "excluded" : "absent";}
		String strata(){
			if(strata.isEmpty()){return ".";}
			final ByteBuilder bb=new ByteBuilder();
			for(String key : strata.keySet()){
				if(bb.length()>0){bb.comma();}
				bb.append(key).colon().append(strata.get(key).intValue());
			}
			return bb.toString();
		}
	}

	private static final String[] FATE_HEADER={
		"schema","stop_id","contig","strand","terminal","longest_start","longest_end","pack_index","model_rank",
		"family_id","representative","identity","score64","tail_length","min_length","min_start_score","raw_stop_score",
		"admitted","rejected_start_score","rejected_length","rejected_limit","soft_rejected","retained","selected","fate",
		"first_pruning_predicates","best_start","best_end","best_raw_start","best_raw_stop","best_inner","best_native",
		"best_bonus","best_total","overlap_basis","overlap_start","overlap_end","max_overlap_same",
		"max_overlap_opposite","observed_selected_overlaps"};
	private static final String[] ROW_HEADER={"schema","input_row","type","contig","strand","start","end","phase",
		"table","partial5","partial3","stratum","attributes"};
	private static final String[] COMPARISON_HEADER={"schema","contig","strand","terminal","reference_raw_rows",
		"reference_primary_rows","baseline_raw_rows","baseline_primary_rows","candidate_raw_rows","candidate_primary_rows",
		"change","reference_rows","baseline_rows","candidate_rows"};
	private static final String STOP_HEADER="genome\tstop_id\tcontig\tstrand\tterminal\treference_raw_rows"
		+"\treference_primary_rows\treference_state\treference_strata\tadmitted\trejected_start_score"
		+"\trejected_length\trejected_limit\tsoft_rejected\tretained\tselected\tfate\tfirst_pruning_predicates"
		+"\tfamily_id\tidentity\tscore64\tlongest_start\tlongest_end\tbest_start\tbest_end\tbest_raw_start"
		+"\tbest_raw_stop\tbest_inner\tbest_native\tbest_bonus\tbest_total\toverlap_basis\toverlap_start"
		+"\toverlap_end\tobserved_selected_overlaps";
	private static final Comparator<Fate> ORDER=new Comparator<Fate>(){
		@Override public int compare(Fate a, Fate b){return a.key.compareTo(b.key);}
	};

	private ProteinSupportFateReport(){}
}
