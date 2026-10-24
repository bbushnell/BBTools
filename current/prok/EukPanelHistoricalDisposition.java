package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Applies the reviewed historical evidence classes and exposes physical-locus links.
 * This conserved annotation ledger is an addendum proposal, never frozen truth.
 * Existing BELOW_GA loci and opposite-strand links remain explicit so a later
 * materialization cannot silently duplicate or promote an existing locus.
 * @author Brian, Raiden */
public final class EukPanelHistoricalDisposition {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelHistoricalDisposition.class, false);
		String truth=null, evidence=null, models=null, history=null, hash=null, out=null, inputs=null;
		int expected=-1, expectedTruth=-1, expectedHistory=-1;
		boolean materialize=false;
		int expectedAssemblies=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("truth")){truth=v;}else if(k.equals("evidence")){evidence=v;}
			else if(k.equals("models")){models=v;}else if(k.equals("history")){history=v;}
			else if(k.equals("truthsha80")){hash=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else if(k.equals("expectedtruth")){expectedTruth=Integer.parseInt(v);}
			else if(k.equals("expectedhistory")){expectedHistory=Integer.parseInt(v);}
			else if(k.equals("materialize")){materialize=parse.Parse.parseBoolean(v);}
			else if(k.equals("inputs")){inputs=v;}else if(k.equals("expectedassemblies")){expectedAssemblies=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(truth!=null && evidence!=null && models!=null && history!=null && out!=null && SeedOffsetTable.validSha80(hash)
			&& expected>0 && expectedTruth>0 && expectedHistory>=expected, "Bound inputs and explicit conservation counts required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		run(truth, evidence, models, history, hash, out, expectedTruth, expectedHistory, expected, materialize, inputs, expectedAssemblies);
		Shared.closeStream(pp.outstream);
	}
	public static void run(String truth, String evidence, String models, String history, String hash,
		String out, int expectedTruth, int expectedHistory, int expected){
		run(truth, evidence, models, history, hash, out, expectedTruth, expectedHistory, expected, false);
	}
	public static void run(String truth, String evidence, String models, String history, String hash,
		String out, int expectedTruth, int expectedHistory, int expected, boolean materialize){
		run(truth, evidence, models, history, hash, out, expectedTruth, expectedHistory, expected, materialize, null, -1);
	}
	public static void run(String truth, String evidence, String models, String history, String hash,
		String out, int expectedTruth, int expectedHistory, int expected, boolean materialize, String inputs, int expectedAssemblies){
		assert(expected>0 && expectedTruth>0):"The disposition ledger must conserve declared baseline and review universes";
		require(!materialize || inputs!=null && expectedAssemblies>0,"Materialization needs the complete panel roster, including assemblies absent from baseline truth");
		final EukPanelHistoricalDisposition x=new EukPanelHistoricalDisposition();
		x.loadModels(models);x.loadTruth(truth, expectedTruth);x.loadEvidence(evidence, hash, expected);
		x.checkHistory(history, expectedHistory);x.resolve();x.write(out);
		if(materialize){new EukPanelHistoricalTruthV2(x, inputs, expectedAssemblies).write(history, out);}
	}
	void loadModels(String path){
		assert(models.isEmpty()):"One sealed roster binds every recommended CM to its family and GA";
		EukPanelHistoricalEvidence.read(path, EukPanelHistoricalReview.MODEL_HEADER, p->{
			require(p.terms()==8, "Model schema changed");final Model m=new Model(p);
			require(m.length>0 && Double.isFinite(m.ga) && models.put(m.id, m)==null, "Invalid or duplicate CM metadata");
		});
	}
	void loadTruth(String path, int expected){
		assert(baseline.isEmpty()):"Existing admitted and BELOW_GA loci both participate in the immutable-baseline join";
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		EukPanelHistoricalEvidence.read(path, EukPanelTruthFreeze.HEADER, p->{
			require(p.terms()==35 && p.termEquals("TRUTH_V1", 32) && p.termEquals(EukPanelTruthDisposition.LEGACY_POLICY, 34), "Only frozen euk truth_v1 is this addendum's baseline");
			final Locus l=new Locus(p);
			require(ids.add(l.id) && (l.admitted || p.termEquals("BELOW_GA", 26)), "Duplicate locus or unknown baseline admission");
			require(l.admitted==EukPanelLocusDiagnostics.bit(p, 14), "Baseline admission must equal its GA evidence");
			ArrayList<Locus> list=baseline.get(l.partition);
			if(list==null){list=new ArrayList<Locus>();baseline.put(l.partition, list);}list.add(l);baselineRows.add(l);
		});
		require(ids.size()==expected, "Baseline locus count changed");
		for(ArrayList<Locus> list:baseline.values()){if(list!=null){Collections.sort(list, LOCUS_ORDER);}}
	}
	void loadEvidence(String path, String hash, int expected){
		assert(rows.isEmpty()):"Each reviewed annotation is consumed exactly once";
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		EukPanelHistoricalEvidence.read(path, INPUT_HEADER, p->{
			require(p.terms()==113 && p.termEquals(hash, 29) && p.termEquals("TRUTH_V2_CANDIDATE_ADDENDUM_NOT_FROZEN", 30), "Evidence schema, baseline binding or status changed");
			final Row r=new Row(p);
			require(ids.add(r.id) && bySource.put(r.sourceKey, r)==null, "Duplicate historical annotation/source identity");
			require(p.termEquals("1", 8) && SeedOffsetTable.validSha80(p.parseString(10)), "Only original eligible, provenance-bound masks were reviewed");
			r.hit=readHit(p, 79, r, false);r.opposite=readHit(p, 95, r, true);
			require(EukPanelLocusDiagnostics.bit(p, 26)==(r.opposite!=null), "Opposite support flag must match its actual qualifying hit");
			if(r.classification.equals("SUPPORTED_GA")){
				require(r.hit!=null && r.hit.ga, "SUPPORTED_GA requires same-strand GA support");r.action="ADMITTED";
			}else if(r.classification.equals("SUPPORTED_SUBGA")){
				require(r.hit!=null && !r.hit.ga && r.hit.subga, "SUPPORTED_SUBGA requires intact significant below-GA evidence");r.action="BELOW_GA";
			}else{
				require(r.classification.equals("UNSUPPORTED") && r.hit==null, "Unknown class or unsupported row with a recommended hit");
				r.action=r.opposite==null?"LIFT_MASK":"ORIENTATION_ERROR_NO_NEW_LOCUS";
			}
			rows.add(r);
		});
		require(rows.size()==expected, "Disposition must conserve every reviewed annotation");
	}
	Hit readHit(LineParser1 p, int offset, Row r, boolean opposite){
		assert(offset==79 || offset==95):"Hit blocks follow the reviewed113-column evidence schema";
		if(p.termEquals("NA", offset)){
			for(int i=offset; i<offset+16; i++){require(p.termEquals("NA", i), "Absent hit must have sixteen NA fields");}return null;
		}
		final Model m=models.get(p.parseString(offset));
		require(m!=null && m.family.equals(r.family), "Recommended CM must belong to the annotation's fold");
		final Hit h=new Hit(p, offset, m);
		require(h.end<=r.contigLength && h.start>=r.windowStart && h.end<=r.windowEnd, "Hit escaped the original genomic window");
		require(h.strand.equals(opposite?other(r.strand):r.strand), "Evidence strand does not match its declared block");
		final long overlap=Math.max(0, Math.min(h.end, r.end)-Math.max(h.start, r.start)+1);
		require(overlap==p.parseLong(offset+13) && EukPanelHistoricalEvidence.reciprocal(h.start, h.end, r.start, r.end), "Qualifying hit must preserve reciprocal-half annotation overlap");
		require(h.ga || h.subga, "Hit does not satisfy GA or full-model untruncated subGA support");return h;
	}
	void checkHistory(String path, int expected){
		assert(!rows.isEmpty()):"Historical inventory validates the original annotation bytes, not just IDs";
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);final int[] counts=new int[2];
		EukPanelHistoricalEvidence.read(path, EukPanelHistoricalMask.OUTPUT_HEADER, p->{
			require(p.terms()==14, "Historical inventory schema changed");counts[0]++;
			final String key=p.parseString(0)+'\t'+p.parseString(8)+'\t'+p.parseLong(10);
			require(seen.add(key), "Historical inventory contains a duplicate source identity");
			final Row r=bySource.get(key);final boolean eligible=EukPanelLocusDiagnostics.bit(p, 7);
			require(eligible==(r!=null), "Every eligible mask, and only eligible masks, must have sensitive review");
			if(r==null){return;}counts[1]++;
			final LineParser1 q=new LineParser1('\t');q.set(r.raw);
			for(int i=0; i<14; i++){
				p.setBounds(i);q.setBounds(i+1);final int n=p.b()-p.a();
				require(n==q.b()-q.a(), "Original mask field length changed for "+r.id);
				for(int j=0; j<n; j++){require(p.line()[p.a()+j]==q.line()[q.a()+j], "Original mask field bytes changed for "+r.id);}
			}
		});
		require(counts[0]==expected && counts[1]==rows.size(), "Historical inventory/review conservation failed");
	}
	void resolve(){
		assert(groups.isEmpty()):"Physical coordinate groups are computed once without changing reference geometry";
		for(Row r:rows){
			if(r.hit!=null){
				r.sameLinks=links(r, r.hit);final String key=key(r.assembly, r.seqid, r.hit.strand, r.family)+'\t'+r.hit.start+'\t'+r.hit.end;
				Group g=byGeometry.get(key);if(g==null){g=new Group("HIST_V2_"+r.id, r);byGeometry.put(key, g);groups.add(g);}
				g.rows.add(r);r.group=g;
			}
			if(r.opposite!=null){r.oppositeLinks=links(r, r.opposite);}
		}
		for(int i=0; i<groups.size(); i++){
			final Group a=groups.get(i);final Row x=a.seed;
			for(int j=i+1; j<groups.size(); j++){
				final Group b=groups.get(j);final Row y=b.seed;
				if(key(x.assembly, x.seqid, x.strand, x.family).equals(key(y.assembly, y.seqid, y.strand, y.family))
					&& EukPanelHistoricalEvidence.reciprocal(x.hit.start, x.hit.end, y.hit.start, y.hit.end)){
					a.overlapping.add(b.id);b.overlapping.add(a.id);
				}
			}
		}
	}
	ArrayList<Locus> links(Row r, Hit h){
		assert(h!=null):"Existing-locus lookup uses the hit's strand, including opposite-strand support";
		final ArrayList<Locus> out=new ArrayList<Locus>(), list=baseline.get(key(r.assembly, r.seqid, h.strand, r.family));
		if(list!=null){for(Locus l:list){if(l.start>h.end){break;}
			if(EukPanelHistoricalEvidence.reciprocal(h.start, h.end, l.start, l.end)){out.add(l);}
		}}return out;
	}
	void write(String prefix){
		assert(!rows.isEmpty()):"Disposition output retains unsupported annotations as explicit decisions";
		RrnaResourceIO.fresh(prefix, ".dispositions.tsv.gz", ".delta.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".dispositions.tsv.gz"), sum=RrnaResourceIO.writer(prefix+".delta.tsv");
		final ByteBuilder b=new ByteBuilder();final ObjectMap<String,long[]> bins=new ObjectMap<String,long[]>(String.class, long[].class);
		try{
			out.println(INPUT_HEADER+"\tdisposition\texact_coordinate_group\texact_group_annotations\tother_reciprocal_groups\texisting_admitted\texisting_below_ga\topposite_existing_admitted\topposite_existing_below_ga\tmaterialization_status");
			for(Row r:rows){
				final long[] n=bins.get(r.family)==null?new long[10]:bins.get(r.family);bins.put(r.family, n);n[0]++;
				if(r.classification.equals("SUPPORTED_GA")){n[1]++;}else if(r.classification.equals("SUPPORTED_SUBGA")){n[2]++;}else{n[3]++;if(r.opposite!=null){n[4]++;}}
				final int admitted=count(r.sameLinks, true), below=count(r.sameLinks, false);
				if(admitted>0){n[5]++;}if(below>0){n[6]++;}
				if(r.group!=null && r==r.group.seed){n[7]++;if(admitted==0 && below==0){n[8]++;}}
				if(r.group!=null && !r.group.overlapping.isEmpty()){n[9]++;}
				b.clear().append(r.raw).tab().append(r.action).tab().append(r.group==null?"NONE":r.group.id).tab().append(r.group==null?0:r.group.rows.size()).tab();
				appendGroups(b, r.group);b.tab();appendLinks(b, r.sameLinks, true);b.tab();appendLinks(b, r.sameLinks, false);
				b.tab();appendLinks(b, r.oppositeLinks, true);b.tab();appendLinks(b, r.oppositeLinks, false);
				out.print(b.tab().append("ADDENDUM_REVIEW_NOT_FROZEN").nl());
			}
			sum.println("fold_family\tannotations\tsupported_ga\tsupported_subga\tunsupported\tunsupported_opposite_only\tannotations_linked_admitted\tannotations_linked_below_ga\texact_supported_coordinate_groups\tgroups_without_baseline_link\tannotations_with_other_reciprocal_groups");
			long total=0;for(String key:EukPanelCmProposal.sorted(bins)){
				final long[] n=bins.get(key);require(n[0]==n[1]+n[2]+n[3] && n[4]<=n[3], "Class disposition must conserve each family");
				b.clear().append(key);for(long v:n){b.tab().append(v);}sum.print(b.nl());total+=n[0];
			}
			require(total==rows.size(), "Published family totals must conserve reviewed annotations");
		}finally{RrnaResourceIO.finish(out, sum);}
		System.err.println("HISTORICAL_DISPOSITION_PASS annotations="+rows.size()+" exact_supported_groups="+groups.size()+" baseline_unchanged=true frozen=false");
	}
	static int count(ArrayList<Locus> list, boolean admitted){int n=0;for(Locus l:list){if(l.admitted==admitted){n++;}}return n;}
	static void appendLinks(ByteBuilder b, ArrayList<Locus> list, boolean admitted){
		int n=0;for(Locus l:list){if(l.admitted==admitted){if(n++>0){b.append(';');}b.append(l.id);}}if(n==0){b.append("NONE");}
	}
	static void appendGroups(ByteBuilder b, Group g){
		if(g==null || g.overlapping.isEmpty()){b.append("NONE");return;}
		for(int i=0; i<g.overlapping.size(); i++){if(i>0){b.append(';');}b.append(g.overlapping.get(i));}
	}
	static String key(String assembly, String seqid, String strand, String family){return assembly+'\t'+seqid+'\t'+strand+'\t'+family;}
	static String other(String strand){require(strand.equals("+") || strand.equals("-"), "Resolved strand required");return strand.equals("+")?"-":"+";}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Model{
		Model(LineParser1 p){id=p.parseString(1);family=p.parseString(2);length=p.parseInt(4);ga=p.parseDouble(5);}
		final String id, family;final int length;final double ga;
	}
	static final class Locus{
		Locus(LineParser1 p){raw=p.line().clone();id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);family=p.parseString(6);bin=p.parseString(7);compartment=p.parseString(8);partition=key(assembly, seqid, strand, family);
			start=p.parseLong(4);end=p.parseLong(5);admitted=p.termEquals("ADMITTED", 26);
			require(start>0 && end>=start && (p.termEquals("+", 3) || p.termEquals("-", 3)), "Baseline geometry/strand invalid");}
		final byte[] raw;final String id, partition, assembly, seqid, strand, family, bin, compartment;final long start, end;final boolean admitted;
	}
	static final class Hit{
		Hit(LineParser1 p, int c, Model m){
			start=p.parseLong(c+2);end=p.parseLong(c+3);strand=p.parseString(c+4);
			// Infernal E-values use scientific notation; the native decimal-only parser rejects exponents.
			final int from=p.parseInt(c+5), to=p.parseInt(c+6);final double bits=p.parseDouble(c+9), evalue=Double.parseDouble(p.parseString(c+10));
			require(p.termEquals("no", c+8) || p.termEquals("5'", c+8) || p.termEquals("3'", c+8) || p.termEquals("5'&3'", c+8), "Unknown Infernal truncation flag");
			require(start>0 && end>=start && p.parseLong(c+1)>0 && from>0 && to>=from && to<=m.length && p.parseInt(c+7)==m.length
				&& p.parseDouble(c+11)==m.ga && Double.isFinite(bits) && bits>=0 && Double.isFinite(evalue) && evalue>=0, "Hit geometry/model/score binding invalid");
			ga=bits>=m.ga;subga=!ga && from==1 && to==m.length && p.termEquals("no", c+8) && evalue<=0.01;
		}
		final long start, end;final String strand;final boolean ga, subga;
	}
	static final class Row{
		Row(LineParser1 p){
			raw=p.line().clone();id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);family=p.parseString(6);
			start=p.parseLong(4);end=p.parseLong(5);contigLength=p.parseLong(15);windowStart=p.parseLong(17);windowEnd=p.parseLong(18);
			sourceKey=assembly+'\t'+p.parseString(9)+'\t'+p.parseLong(11);classification=p.parseString(23);
			require(start>0 && end>=start && windowStart>0 && windowStart<=start && windowEnd>=end && windowEnd<=contigLength
				&& p.parseLong(20)==end-start+1 && p.parseLong(11)>0 && (strand.equals("+") || strand.equals("-")), "Original annotation/window geometry invalid");
		}
		final byte[] raw;final String id, assembly, seqid, strand, family, sourceKey, classification;
		final long start, end, contigLength, windowStart, windowEnd;
		String action;Hit hit, opposite;Group group;
		ArrayList<Locus> sameLinks=new ArrayList<Locus>(), oppositeLinks=new ArrayList<Locus>();
	}
	static final class Group{
		Group(String id_, Row r){id=id_;seed=r;}
		final String id;final Row seed;final ArrayList<Row> rows=new ArrayList<Row>();final ArrayList<String> overlapping=new ArrayList<String>();
	}
	static final Comparator<Locus> LOCUS_ORDER=(a,b)->{final int c=Long.compare(a.start, b.start);return c==0?a.id.compareTo(b.id):c;};
	public static final String INPUT_HEADER=EukPanelHistoricalReview.WINDOW_HEADER+"\tannotation_length\tproduct\tannotation_group\trecommended_class\tproposed_action\treason\topposite_strand_flag\texisting_loci\traw_table\tbaseline_truth_sha80\tcandidate_status"
		+EukPanelHistoricalEvidence.columns("same")+EukPanelHistoricalEvidence.columns("subga")+EukPanelHistoricalEvidence.columns("opposite")
		+EukPanelHistoricalEvidence.columns("recommended")+EukPanelHistoricalEvidence.columns("opposite_support")+"\tsame_max_overlap_ties\topposite_max_overlap_ties";
	final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class, Model.class);
	@SuppressWarnings("unchecked")
	final ObjectMap<String,ArrayList<Locus>> baseline=new ObjectMap(String.class, ArrayList.class);
	final ObjectMap<String,Row> bySource=new ObjectMap<String,Row>(String.class, Row.class);
	final ObjectMap<String,Group> byGeometry=new ObjectMap<String,Group>(String.class, Group.class);
	final ArrayList<Row> rows=new ArrayList<Row>();final ArrayList<Group> groups=new ArrayList<Group>();
	final ArrayList<Locus> baselineRows=new ArrayList<Locus>();
}
