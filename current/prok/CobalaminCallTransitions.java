package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import shared.Shared;
import structures.ByteBuilder;

/** Compares preserved union assignments, exact FP survival and spatial gene-row
 * changes. Overlap-only links never become exact survival or causal claims.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCallTransitions {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminCallTransitions.class, false);
		String before=null, after=null, panel=null, contigs=null, changes=null, out=null;
		int beforeCount=-1, afterCount=-1;
		for(String arg:pp.args){
			final int equals=arg.indexOf('=');require(equals>0, "Use flag=value");
			final String key=arg.substring(0, equals).toLowerCase(Locale.ROOT), value=arg.substring(equals+1);
			if(key.equals("before")){before=value;}else if(key.equals("after")){after=value;}
			else if(key.equals("panel")){panel=value;}else if(key.equals("contigs")){contigs=value;}
			else if(key.equals("changes")){changes=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("beforecount")){beforeCount=Integer.parseInt(value);}else if(key.equals("aftercount")){afterCount=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(before!=null && after!=null && panel!=null && contigs!=null && out!=null && beforeCount>=0 && afterCount>=0,
			"Bind both native call tables, the same frozen panel, ownership and expected call counts");
		Shared.setThreads(1);
		try{
			final CobalaminCallTransitions x=new CobalaminCallTransitions();x.panel(panel, contigs);
			final Run a=x.calls(before, beforeCount), z=x.calls(after, afterCount);
			x.write(a, z, changes, out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private void panel(String path, String contigs){
		read(contigs, CobalaminPanelPrep.CONTIG_HEADER, p->{
			final String seq=p.parseString(2), genome=p.parseString(0);final int length=p.parseInt(3);
			require(!seq.isEmpty() && !genome.isEmpty() && length>0 && !owners.contains(seq), "Contig ownership must be unique and explicit");
			owners.put(seq, genome);lengths.put(seq, length);
		});
		read(path, CobalaminPanelPrep.HEADER, p->{
			final String genome=p.parseString(1), seq=p.parseString(2), key=genome+'\t'+p.parseString(0);
			final String source=p.parseString(3), member=p.parseString(4);
			require(genome.equals(owners.get(seq)) && known(source) && !truth.contains(key), "Truth identity must retain its original genome and source family");
			final boolean full=!member.equals(CobalaminPanelPrep.UNRESOLVED);
			require(!full || known(member) || member.equals(CobalaminPanelPrep.CLAN_ONLY), "Unknown reviewed truth membership");
			final int start=p.parseInt(7), stop=p.parseInt(8);final String strand=p.parseString(9);
			require((strand.equals("+") || strand.equals("-")) && (full ? start>0 && stop>=start && stop<=lengths.get(seq) : start==0 && stop==0),
				"Admitted truth must retain its exact full-core geometry");
			truth.put(key, new Truth(source, member, full ? new Locus(seq, "B12_union", start, stop, strand.charAt(0)) : null));
		});
		require(owners.size()>0 && truth.size()>0, "An empty panel cannot validate saved assignments");
	}
	private Run calls(String path, int expected){
		final Run run=new Run();
		read(path, NcrnaReviewUnion.PANEL_CALL_HEADER, p->{
			if(!p.termEquals("FULL_CM", 0) || !p.termEquals("union", 1)){return;}
			final Call c=new Call(p);require(c.genome.equals(owners.get(c.seqid)) && c.stop<=lengths.get(c.seqid), "Call is outside its original source contig");
			require(!run.ids.contains(c.id) && !run.exact.contains(c.key()), "Duplicate saved calls require explicit event accounting");
			run.ids.put(c.id, c);run.exact.put(c.key(), c);run.calls.add(c);
			Group group=run.contigs.get(c.seqid);if(group==null){group=new Group();run.contigs.put(c.seqid, group);}group.calls.add(c);
			if(c.outcome.equals("TP")){
				final String key=c.genome+'\t'+c.target;final Truth target=truth.get(key);
				require(target!=null && target.locus!=null && !run.assigned.contains(key)
					&& NcrnaGffFamilyGrader.overlapsReciprocal(c.locus, target.locus, .5f),
					"Saved TPs must retain one-to-one, same-strand reciprocal-half assignments to admitted truth");
				c.source=target.source;c.member=target.member;run.assigned.put(key, c);run.tp++;
			}else{
				require(c.target.equals("NA"), "An unassigned call cannot carry a native truth assignment");
				if(c.outcome.equals("APPARENT_FP")){run.fp++;}else{run.unscored++;}
			}
		});
		require(run.calls.size()==expected && expected==run.tp+run.fp+run.unscored, "Saved union call denominator changed: "+path);
		return run;
	}
	private void write(Run before, Run after, String changes, String out){
		RrnaResourceIO.fresh(out, ".summary.tsv", ".fp.tsv", ".sibling_tp.tsv", ".gene_links.tsv", ".tp.tsv", ".labels.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv"), fp=RrnaResourceIO.writer(out+".fp.tsv");
		final ByteStreamWriter tp=RrnaResourceIO.writer(out+".sibling_tp.tsv");final ByteBuilder b=new ByteBuilder();
		int exact=0, overlap=0, absent=0, sibling=0, recovered=0, retained=0, lost=0;
		try{
			fp.println("before_call\tgenome\tseqid\tstart1\tstop1\tstrand\tbefore_family\tstatus\tafter_call\tafter_family\tafter_outcome\tafter_start1\tafter_stop1\toverlap_bp");
			tp.println("before_call\tgenome\ttruth_id\tsource_family\ttruth_member_family\tbefore_family\tstatus\tafter_call\tafter_family\texact_same_span");
			for(Call c:before.calls){
				if(c.outcome.equals("APPARENT_FP")){
					final Call match=after.exact.get(c.key());
					if(match!=null){exact++;fpRow(fp, b, c, match, "EXACT_PRESENT");}
					else{
						int links=0;final Group group=after.contigs.get(c.seqid);
						if(group!=null){for(Call z:group.calls){if(c.strand.equals(z.strand) && overlap(c.start, c.stop, z.start, z.stop)>0){
							links++;fpRow(fp, b, c, z, "OVERLAP_ONLY");
						}}}
						if(links>0){overlap++;}else{absent++;fpRow(fp, b, c, null, "ABSENT");}
					}
				}
				if(c.outcome.equals("TP") && !c.family.equals("RF00174") && c.source.equals("RF00174")){
					sibling++;final Call match=after.assigned.get(c.genome+'\t'+c.target);
					final String state;
					if(match==null){lost++;state="LOST_NATIVE_ASSIGNMENT";}
					else if(match.family.equals("RF00174")){recovered++;state="RECOVERED_BY_RF00174";}
					else{retained++;state="RETAINED_BY_SIBLING";}
					tp.print(b.clear().append(c.id).tab().append(c.genome).tab().append(c.target).tab().append(c.source).tab().append(c.member).tab().append(c.family)
						.tab().append(state).tab().append(match==null?"NA":match.id).tab().append(match==null?"NA":match.family)
						.tab().append(match==null?"NA":Boolean.toString(c.key().equals(match.key()))).nl());
				}
			}
			require(exact+overlap+absent==before.fp && recovered+retained+lost==sibling,
				"Transition bins must conserve original FP calls and RF00174-assigned sibling TPs");
			summary.println("metric\tcount");
			final String[] labels={"before_calls", "before_TP", "before_FP", "before_unscored", "after_calls", "after_TP", "after_FP", "after_unscored",
				"before_FP_exact_present", "before_FP_overlap_only", "before_FP_absent", "before_sibling_TP_on_RF00174_truth",
				"recovered_by_RF00174", "retained_by_sibling", "lost_native_assignment"};
			final int[] values={before.calls.size(), before.tp, before.fp, before.unscored, after.calls.size(), after.tp, after.fp, after.unscored,
				exact, overlap, absent, sibling, recovered, retained, lost};
			for(int i=0; i<labels.length; i++){summary.print(b.clear().append(labels[i]).tab().append(values[i]).nl());}
		}finally{RrnaResourceIO.finish(summary, fp, tp);}
		allTransitions(before, after, out);
		if(changes!=null){geneLinks(changes, after, out+".gene_links.tsv");}
		System.err.println("COBALAMIN_CALL_TRANSITIONS_PASS before="+before.calls.size()+" after="+after.calls.size()+"; exact survival and native assignments, not causal inference");
	}
	/** Attribute losses to the original emitting label while retaining the same
	 * one-to-one truth identity. A changed label is not a lost biological locus. */
	private static void allTransitions(Run before, Run after, String prefix){
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".tp.tsv"), labels=RrnaResourceIO.writer(prefix+".labels.tsv");
		final ByteBuilder b=new ByteBuilder();final int[][] counts=new int[3][14];int accounted=0;
		try{
			detail.println("before_call\tgenome\ttruth_id\tsource_family\ttruth_member_family\tbefore_family\tstatus\tafter_call\tafter_family\texact_same_span");
			for(Call c:after.calls){counts[familyIndex(c.family)][3+outcomeIndex(c.outcome)]++;}
			for(Call c:before.calls){
				final int[] n=counts[familyIndex(c.family)];n[outcomeIndex(c.outcome)]++;
				if(c.outcome.equals("TP")){
					accounted++;final Call z=after.assigned.get(c.genome+'\t'+c.target);final String state;
					if(z==null){n[8]++;state="LOST_NATIVE_ASSIGNMENT";}
					else if(c.family.equals(z.family)){n[6]++;state="RETAINED_SAME_LABEL";}
					else{n[7]++;state="RETAINED_CHANGED_LABEL";}
					if(!c.family.equals("RF00174") && c.source.equals("RF00174")){
						if(z==null){n[10]++;}else if(z.family.equals("RF00174")){n[9]++;}
					}
					detail.print(b.clear().append(c.id).tab().append(c.genome).tab().append(c.target).tab().append(c.source).tab().append(c.member).tab().append(c.family)
						.tab().append(state).tab().append(z==null?"NA":z.id).tab().append(z==null?"NA":z.family)
						.tab().append(z==null?"NA":Boolean.toString(c.key().equals(z.key()))).nl());
				}else if(c.outcome.equals("APPARENT_FP")){
					if(after.exact.contains(c.key())){n[11]++;}else{
						final Group group=after.contigs.get(c.seqid);boolean overlaps=false;
						if(group!=null){for(Call z:group.calls){if(c.strand.equals(z.strand) && overlap(c.start, c.stop, z.start, z.stop)>0){overlaps=true;break;}}}
						n[overlaps?12:13]++;
					}
				}
			}
			require(accounted==before.tp, "Every original union TP must retain a truth-identity transition");
			labels.println("family\tbefore_TP\tbefore_FP\tbefore_unscored\tafter_TP\tafter_FP\tafter_unscored\tTP_same_label\tTP_changed_label\tTP_lost_assignment\tsibling_RF00174_truth_recovered_by_RF00174\tsibling_RF00174_truth_lost\tFP_exact_present\tFP_overlap_only\tFP_absent");
			for(int i=0; i<counts.length; i++){
				final int[] n=counts[i];require(n[0]==n[6]+n[7]+n[8] && n[1]==n[11]+n[12]+n[13], "Emitting-label bins must conserve original TPs and FPs");
				b.clear().append(FAMILIES[i]);for(int value:n){b.tab().append(value);}labels.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(detail, labels);}
	}
	private static int familyIndex(String family){for(int i=0; i<FAMILIES.length; i++){if(FAMILIES[i].equals(family)){return i;}}throw new IllegalArgumentException("Unknown emitted family: "+family);}
	private static int outcomeIndex(String outcome){return outcome.equals("TP")?0:outcome.equals("APPARENT_FP")?1:2;}
	private static void fpRow(ByteStreamWriter out, ByteBuilder b, Call c, Call z, String state){
		out.print(b.clear().append(c.id).tab().append(c.genome).tab().append(c.seqid).tab().append(c.start).tab().append(c.stop).tab().append(c.strand)
			.tab().append(c.family).tab().append(state).tab().append(z==null?"NA":z.id).tab().append(z==null?"NA":z.family)
			.tab().append(z==null?"NA":z.outcome).tab().append(z==null?"NA":Integer.toString(z.start)).tab().append(z==null?"NA":Integer.toString(z.stop))
			.tab().append(z==null?0:overlap(c.start, c.stop, z.start, z.stop)).nl());
	}
	/** All removed rows, including changed boundaries/scores, remain visible. The
	 * native delta retains added rows separately; these links do not classify loss. */
	private void geneLinks(String path, Run after, String output){
		final ByteStreamWriter out=RrnaResourceIO.writer(output);final ByteBuilder b=new ByteBuilder();
		final GffCallSetDelta classifier=new GffCallSetDelta();final long[] line={1};
		try{
			out.println("delta_line\tchanged_family\tgenome\tb12_call\tb12_family\tb12_outcome\tstrand_relation\toverlap_bp\t"+GffCallSetDelta.CHANGE_HEADER);
			read(path, GffCallSetDelta.CHANGE_HEADER, p->{
				line[0]++;require(p.termEquals("REMOVED", 0) || p.termEquals("ADDED", 0), "Unknown native delta direction");
				if(!p.termEquals("REMOVED", 0)){return;}
				p.setBounds(1);final byte[] gff=Arrays.copyOfRange(p.line(), p.a(), p.line().length);
				final String family=classifier.family(gff);if(NcrnaFamilyRegistry.isCobalamin(family)){return;}
				final String seq=p.parseString(1), genome=owners.get(seq), strand=p.parseString(7);
				final int start=p.parseInt(4), stop=p.parseInt(5);
				require(genome!=null && start>0 && stop>=start && stop<=lengths.get(seq) && (strand.equals("+") || strand.equals("-")),
					"Removed gene row must bind the same frozen panel geometry");
				final Group group=after.contigs.get(seq);int links=0;
				if(group!=null){for(Call c:group.calls){final int bp=overlap(start, stop, c.start, c.stop);if(bp<=0){continue;}links++;
					out.print(b.clear().append(line[0]).tab().append(family).tab().append(genome).tab().append(c.id).tab().append(c.family).tab().append(c.outcome)
						.tab().append(strand.equals(c.strand)?"SAME":"OPPOSITE").tab().append(bp).tab().append(p.line()).nl());
				}}
				if(links==0){out.print(b.clear().append(line[0]).tab().append(family).tab().append(genome).append("\tNA\tNA\tNA\tNO_B12_OVERLAP\t0\t").append(p.line()).nl());}
			});
		}finally{RrnaResourceIO.finish(out);}
	}
	private static void read(String path, String header, java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		p.set(header.getBytes(java.nio.charset.StandardCharsets.US_ASCII));final int columns=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==columns, "Native table schema changed: "+path);action.accept(p);
		}}finally{require(!in.close(), "Native table read failed: "+path);}
	}
	private static boolean known(String family){return family.equals("RF00174") || family.equals("RF01482") || family.equals("RF01689");}
	private static int overlap(int a, int z, int b, int y){return Math.max(0, Math.min(z, y)-Math.max(a, b)+1);}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Call {
		Call(LineParser1 p){
			id=p.parseString(2);genome=p.parseString(3);seqid=p.parseString(4);family=p.parseString(5);start=p.parseInt(6);stop=p.parseInt(7);
			strand=p.parseString(8);model=p.parseString(9);outcome=p.parseString(10);target=p.parseString(11);
			require(!id.isEmpty() && known(family) && start>0 && stop>=start && (strand.equals("+") || strand.equals("-"))
				&& model.startsWith("cobalamin_"+family+'_')
				&& (outcome.equals("TP") || outcome.equals("APPARENT_FP") || outcome.equals("UNSCORED")), "Malformed or unattributed saved union call");
			locus=new Locus(seqid, "B12_union", start, stop, strand.charAt(0));
		}
		String key(){return seqid+'\t'+strand+'\t'+start+'\t'+stop;}
		final String id, genome, seqid, family, strand, model, outcome, target;final int start, stop;final Locus locus;String source, member;
	}
	private static final class Truth {Truth(String s, String m, Locus l){source=s;member=m;locus=l;}final String source, member;final Locus locus;}
	private static final class Group {final ArrayList<Call> calls=new ArrayList<Call>();}
	private static final class Run {
		final ArrayList<Call> calls=new ArrayList<Call>();
		final ObjectMap<String,Call> ids=new ObjectMap<String,Call>(String.class, Call.class), exact=new ObjectMap<String,Call>(String.class, Call.class), assigned=new ObjectMap<String,Call>(String.class, Call.class);
		final ObjectMap<String,Group> contigs=new ObjectMap<String,Group>(String.class, Group.class);int tp, fp, unscored;
	}
	private final ObjectMap<String,String> owners=new ObjectMap<String,String>(String.class, String.class);
	private final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
	private final ObjectMap<String,Truth> truth=new ObjectMap<String,Truth>(String.class, Truth.class);
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
}
