package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import shared.Shared;
import structures.ByteBuilder;

/** Audits saved union assignments and physical call overlap without changing grades.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelAttribution {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminPanelAttribution.class, false);
		String panel=null, contigs=null, runs=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("panel")){panel=v;}else if(k.equals("contigs")){contigs=v;}else if(k.equals("runs")){runs=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedruns")){expected=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(panel!=null && contigs!=null && runs!=null && out!=null && expected>0, "Bind the accepted panel, contig ownership and complete saved run manifest");
		Shared.setThreads(1);
		try{final CobalaminPanelAttribution x=new CobalaminPanelAttribution();x.contigs(contigs);x.truth(panel);x.runs(runs, expected);x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void contigs(String path){
		read(path, CobalaminPanelPrep.CONTIG_HEADER, p->{
			final String genome=p.parseString(0), seqid=p.parseString(2);final int length=p.parseInt(3);
			require(!genome.isEmpty() && length>0 && !owners.contains(seqid), "A contig must retain one nonempty owner and positive length");
			owners.put(seqid, new Contig(genome, length));
		});
		require(owners.size()>0, "Attribution requires the accepted panel's contig map");
	}
	private void truth(String path){
		read(path, CobalaminPanelPrep.HEADER, p->{
			final String id=p.parseString(0), genome=p.parseString(1), seqid=p.parseString(2), source=p.parseString(3), member=p.parseString(4);
			require(family(source), "Unexpected source subtype in accepted panel");
			check(genome, seqid, p.parseInt(5), p.parseInt(6), p.parseString(9));
			final String key=genome+'\t'+id;
			require(sourceIds.add(key), "Ambiguous source ID needs explicit disambiguation before joining saved assignments: "+key);
			if(member.equals(CobalaminPanelPrep.UNRESOLVED)){
				require(p.parseInt(7)==0 && p.parseInt(8)==0, "Unresolved evidence must not acquire full-core coordinates");return;
			}
			require(member.equals(source) || member.equals(CobalaminPanelPrep.CLAN_ONLY), "Reviewed full cores retain own-subtype or clan-only evidence");
			final Locus l=check(genome, seqid, p.parseInt(7), p.parseInt(8), p.parseString(9));
			final Truth t=new Truth(id, source, member, l);truthById.put(key, t);
			Group g=truthByContig.get(seqid);if(g==null){g=new Group();truthByContig.put(seqid, g);}g.truth.add(t);
		});
	}
	private void runs(String path, int expected){
		final ObjectSet<String> names=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		read(path, RUN_HEADER, p->{
			final String name=p.parseString(0), file=p.parseString(4);final int k=p.parseInt(1), top=p.parseInt(2);final double identity=p.parseDouble(3);
			require(!name.isEmpty() && names.add(name) && paths.add(file) && k>0 && top>0 && identity>0 && identity<=1, "Run names/paths must be unique with explicit settings");
			final Run r=new Run(name, k, top, identity);calls(r, file);runs.add(r);
			if(top==100 && Double.compare(identity, .64)==0 && k>=7 && k<=9){require(selected[k-7]==null, "One top100/.64 run per map K is required");selected[k-7]=r;}
		});
		require(runs.size()==expected && selected[0]!=null && selected[1]!=null && selected[2]!=null, "Consume all runs and the three fixed K7/K8/K9 top100/.64 comparisons");
	}
	private void calls(Run r, String path){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		read(path, NcrnaReviewUnion.PANEL_CALL_HEADER, p->{
			if(!p.termEquals("FULL_CM", 0) || !p.termEquals("union", 1)){return;}
			final String id=p.parseString(2), genome=p.parseString(3), seqid=p.parseString(4), emitted=p.parseString(5), outcome=p.parseString(10), target=p.parseString(11);
			require(ids.add(id) && family(emitted) && (outcome.equals("TP") || outcome.equals("APPARENT_FP") || outcome.equals("UNSCORED")), "Saved union audit requires unique call IDs and known labels/outcomes");
			final Locus l=check(genome, seqid, p.parseInt(6), p.parseInt(7), p.parseString(8));
			final Truth assigned=outcome.equals("TP")?truthById.get(genome+'\t'+target):null;
			require(outcome.equals("TP") ? assigned!=null && NcrnaGffFamilyGrader.overlapsReciprocal(l, assigned.locus, .5f) : target.equals("NA"),
				"Saved union assignment must point to the same strand and reciprocal-half full core, or explicitly no target: "+r.name+" "+id);
			final Call c=new Call(id, genome, emitted, outcome, l, assigned);final Group g=truthByContig.get(seqid);
			if(g!=null){for(Truth t:g.truth){if(NcrnaGffFamilyGrader.overlapsReciprocal(l, t.locus, .5f)){
				if(t.source.equals("RF00174")){c.rf1++;}else{c.other++;}
			}}}
			final String key=key(l);require(!r.physical.contains(key), "Exact duplicate calls require separate event accounting before set comparison: "+r.name+" "+key);
			r.physical.put(key, c);r.calls.add(c);r.loci.add(l);
		});
	}
	private Locus check(String genome, String seqid, int start, int stop, String strand){
		final Contig c=owners.get(seqid);
		require(c!=null && c.genome.equals(genome) && start>0 && stop>=start && stop<=c.length && (strand.equals("+") || strand.equals("-")), "Audit coordinates must remain on their exact source genome/contig/strand");
		return new Locus(seqid, "B12_union", start, stop, strand.charAt(0));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String out){
		RrnaResourceIO.fresh(out, ".summary.tsv", ".calls.tsv", ".overlap.tsv", ".pairs.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv"), detail=RrnaResourceIO.writer(out+".calls.tsv");
		final ByteStreamWriter overlap=RrnaResourceIO.writer(out+".overlap.tsv"), pairs=RrnaResourceIO.writer(out+".pairs.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			summary.println("run\tmap_k\ttop_n\tid_pass\tcalls\tTP\tapparent_FP\tunscored\tsibling_calls\tsibling_overlap_RF00174_full\tsibling_overlap_any_full\tsibling_no_full_overlap\tdirect_cross_subtype_TP\tsibling_assigned_RF00174_TP");
			detail.println("run\tcall_id\tgenome\tseqid\tstart1\tstop1\tstrand\temitted_family\tunion_outcome\tassigned_truth_id\tassigned_source_family\tassigned_member_family\tRF00174_full_overlaps\tother_full_overlaps\tdirect_cross_subtype_TP");
			for(Run r:runs){
				int tp=0, fp=0, unscored=0, sibling=0, rf1=0, any=0, none=0, cross=0, assignedRf1=0;
				for(Call c:r.calls){
					final Truth t=c.assigned;final boolean isSibling=!c.family.equals("RF00174");
					final boolean swapped=t!=null && !t.member.equals(CobalaminPanelPrep.CLAN_ONLY) && !t.member.equals(c.family);
					if(t!=null){tp++;}else if(c.outcome.equals("APPARENT_FP")){fp++;}else{unscored++;}
					if(swapped){cross++;}
					if(isSibling){sibling++;if(c.rf1>0){rf1++;}if(c.rf1+c.other>0){any++;}else{none++;}if(t!=null && t.source.equals("RF00174")){assignedRf1++;}}
					detail.print(b.clear().append(r.name).tab().append(c.id).tab().append(c.genome).tab().append(c.locus.seqid).tab().append(c.locus.start).tab().append(c.locus.stop).tab().append(c.locus.strand)
						.tab().append(c.family).tab().append(c.outcome).tab().append(t==null?"NA":t.id).tab().append(t==null?"NA":t.source).tab().append(t==null?"NA":t.member)
						.tab().append(c.rf1).tab().append(c.other).tab().append(swapped).nl());
				}
				require(tp+fp+unscored==r.calls.size() && any+none==sibling && rf1<=any && cross<=tp, "Attribution conserves every saved call; overlap counts are not extra TPs");
				summary.print(b.clear().append(r.name).tab().append(r.k).tab().append(r.top).tab().append(r.identity, 4).tab().append(r.calls.size()).tab().append(tp).tab().append(fp).tab().append(unscored)
					.tab().append(sibling).tab().append(rf1).tab().append(any).tab().append(none).tab().append(cross).tab().append(assignedRf1).nl());
			}
			overlap.println("run_a\trun_b\tcalls_a\tcalls_b\texact_coordinate_shared\texact_shared_different_subtype\treciprocal_half_one_to_one\tunmatched_a\tunmatched_b");
			pairs.println("run_a\trun_b\tcall_a\tcall_b\texact_coordinates\tfamily_a\tfamily_b");
			for(int i=0; i<3; i++){for(int j=i+1; j<3; j++){
				final Run a=selected[i], z=selected[j];int exact=0, labels=0;
				for(Call c:a.calls){final Call d=z.physical.get(key(c.locus));if(d!=null){exact++;if(!c.family.equals(d.family)){labels++;}}}
				final MatchResult m=NcrnaGffFamilyGrader.matchOneToOne(a.loci, z.loci, .5f);
				require(m.matched>=exact, "Maximum-cardinality overlap includes every exact pair as a feasible matching");
				overlap.print(b.clear().append(a.name).tab().append(z.name).tab().append(a.calls.size()).tab().append(z.calls.size()).tab().append(exact).tab().append(labels)
					.tab().append(m.matched).tab().append(a.calls.size()-m.matched).tab().append(z.calls.size()-m.matched).nl());
				for(int q=0; q<a.calls.size(); q++){final int t=m.callMatchedTruthIdx[q];if(t<0){continue;}
					final Call c=a.calls.get(q), d=z.calls.get(t);
					pairs.print(b.clear().append(a.name).tab().append(z.name).tab().append(c.id).tab().append(d.id).tab().append(key(c.locus).equals(key(d.locus))).tab().append(c.family).tab().append(d.family).nl());
				}
			}}
		}finally{RrnaResourceIO.finish(summary, detail, overlap, pairs);}
		System.err.println("COBALAMIN_PANEL_ATTRIBUTION_PASS runs="+runs.size()+"; saved assignments and interval overlaps, not causal competition evidence");
	}
	private static void read(String path, String header, java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.US_ASCII));final int columns=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==columns, "Exact native schema required: "+path);action.accept(p);}}
		finally{require(!in.close(), "Attribution input read failed: "+path);}
	}
	private static boolean family(String name){for(String f:CobalaminPanelPrep.FAMILIES){if(f.equals(name)){return true;}}return false;}
	private static String key(Locus l){return l.seqid+'\t'+l.strand+'\t'+l.start+'\t'+l.stop;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Contig {Contig(String g, int n){genome=g;length=n;}final String genome;final int length;}
	private static final class Truth {Truth(String i, String s, String m, Locus l){id=i;source=s;member=m;locus=l;}final String id, source, member;final Locus locus;}
	private static final class Group {final ArrayList<Truth> truth=new ArrayList<Truth>();}
	private static final class Call {
		Call(String i, String g, String f, String o, Locus l, Truth t){id=i;genome=g;family=f;outcome=o;locus=l;assigned=t;}
		final String id, genome, family, outcome;final Locus locus;final Truth assigned;int rf1, other;
	}
	private static final class Run {
		Run(String n, int k_, int t, double d){name=n;k=k_;top=t;identity=d;}
		final String name;final int k, top;final double identity;final ArrayList<Call> calls=new ArrayList<Call>();final ArrayList<Locus> loci=new ArrayList<Locus>();
		final ObjectMap<String,Call> physical=new ObjectMap<String,Call>(String.class, Call.class);
	}
	private final ObjectMap<String,Contig> owners=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ObjectMap<String,Truth> truthById=new ObjectMap<String,Truth>(String.class, Truth.class);
	private final ObjectMap<String,Group> truthByContig=new ObjectMap<String,Group>(String.class, Group.class);
	private final ObjectSet<String> sourceIds=new ObjectSet<String>(String.class);
	private final ArrayList<Run> runs=new ArrayList<Run>();private final Run[] selected=new Run[3];
	public static final String RUN_HEADER="run\tmap_k\ttop_n\tid_pass\tcalls";
}
