package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

import map.IntHashMap;
import map.IntHashSet;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Exact selected-path multiset comparison and cobalamin reported-overlap accounting.
 * Adapted from template.A_SampleByteFile and NcrnaBoundaryInvocationSummary.
 * Native byte parsing interns distinct text only; one retained Row per selected
 * result supports order-independent comparison across concurrent caller workers.
 * No genome, neural score, alignment or reference matching is recomputed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSelectionAudit{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSelectionAudit.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			int expected=-1;
			String movements=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("expectedinputs")){expected=Parse.parseIntKMG(b);}
				else if(a.equals("movementout")){movements=b;}
				else{require(parser.parse(arg, a, b), "Unknown selection-audit argument: "+arg);}
			}
			require(parser.out2==null && !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Selection auditing requires two complete traces and fresh output");
			Shared.setThreads(1);
			audit(parser.in1, parser.in2, parser.out1, movements, expected);
			pp.outstream.println("NCRNA_SELECTION_AUDIT_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Writes the small result before failing on a genuine selected-path difference. */
	public static void audit(final String before, final String after, final String out, final int expected){
		audit(before, after, out, null, expected);
	}

	/** Optional physical-call movement counts never change the historical summary schema. */
	public static void audit(final String before, final String after, final String out, final String movements, final int expected){
		require(before!=null && after!=null && out!=null && expected>0 && Tools.testInputFiles(false, true, before, after)
			&& Tools.testOutputFiles(false, false, false, out, movements) && Tools.testForDuplicateFiles(true, before, after, out, movements),
			"Supply in1=OFF_geometry in2=ON_geometry out=fresh_summary expectedinputs=N");
		final NcrnaSelectionAudit tool=new NcrnaSelectionAudit(movements!=null);
		final Data off=tool.read(before), on=tool.read(after);
		require(off.inputs.size()==expected && on.inputs.size()==expected, "Both traces must declare the exact input population");
		for(int ordinal=0; ordinal<expected; ordinal++){
			require(off.inputs.get(ordinal)>=0 && off.inputs.get(ordinal)==on.inputs.get(ordinal),
				"Input ordinal/path differs or is missing at ordinal "+ordinal);
		}
		require(off.contigs.size()==on.contigs.size(), "Complete selected traces must retain the same contig population");
		for(final long key:off.contigs.toArray()){
			require(on.contigs.contains(key), "Selected traces differ in input/read population");
		}
		Collections.sort(off.rows);
		Collections.sort(on.rows);
		final int common=Math.min(off.rows.size(), on.rows.size());
		int first=-1;
		for(int i=0; i<common; i++){
			if(off.rows.get(i).compareTo(on.rows.get(i))!=0){first=i;break;}
		}
		if(first<0 && off.rows.size()!=on.rows.size()){first=common;}
		final boolean equal=first<0 && off.rows.size()==on.rows.size();
		try(Output output=new Output(out)){
			output.writer.println("metric\tvalue");
			final ByteBuilder b=new ByteBuilder();
			metric(output, b, "exact_path_multiset", equal ? 1 : 0);
			metric(output, b, "before_selected_features", off.rows.size());
			metric(output, b, "after_selected_features", on.rows.size());
			metric(output, b, "first_sorted_difference", first);
			writeStats(output, b, "before", off);
			writeStats(output, b, "after", on);
		}
		require(equal, "Selected path identities, scores or multiplicities changed; preserve traces and stop the experiment. Summary: "+out);
		if(movements!=null){
			try(Output output=new Output(movements)){
				output.writer.println("view\tmember\tselected_calls\tmoved5\tmoved3\tmoved_both\tmoved_either");
				writeMovements(output, "before", off);
				writeMovements(output, "after", on);
			}
		}
	}

	private NcrnaSelectionAudit(final boolean movements_){movements=movements_;}

	/** The producer prints each P/Q/T contig block atomically; blocks may arrive in any order. */
	private Data read(final String path){
		final Data data=new Data();
		final ArrayList<Row> active=new ArrayList<Row>();
		final IntHashSet[] affected={new IntHashSet(), new IntHashSet()};
		final LongHashSet pairs=new LongHashSet();
		int ordinal=-1, contig=-1;
		long read=-1;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "#ncrna_selected_geometry\t1", path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line);
				final int source=p.parseInt(1);
				require(source>=0, "Selected geometry requires nonnegative input ordinals");
				if(p.termEquals("I", 0)){
					require(p.terms()==3 && active.isEmpty() && data.inputs.get(source)<0, "Input paths must be unique declarations between complete contigs");
					data.inputs.put(source, symbol(2));
					continue;
				}
				require(data.inputs.get(source)>=0, "A contig must follow its exact input declaration");
				final long readID=p.parseLong(2);
				final int name=symbol(3);
				require(readID>=0 && readID<=0xffffffffL, "Selection audit supports unsigned32-bit per-input read IDs");
				if(active.isEmpty()){
					ordinal=source;
					read=readID;
					contig=name;
				}
				require(source==ordinal && readID==read && name==contig, "Contig blocks must not interleave before their conservation row");
				if(p.termEquals("P", 0)){
					require(p.terms()==15 && p.parseInt(4)==active.size(), "Selected local indices must be dense and unique within the contig");
					final Row row=new Row(source, readID, name, p.parseInt(5), p.parseInt(6), p.parseInt(7),
						symbol(8), symbol(9), p.parseInt(10), p.parseInt(11), p.parseInt(12), p.parseInt(13), p.parseInt(14));
					require((row.strand==0 || row.strand==1) && row.frame>=0 && row.frame<3 && row.start>=0
						&& row.stop>=row.start && row.reportStart>=0 && row.reportStop>=row.reportStart,
						"Selected rows require physical ordered intervals, strand and frame");
					active.add(row);
					data.rows.add(row);
					if(row.start!=row.reportStart || row.stop!=row.reportStop){data.moved++;}
					if(movements){countMovement(data, row);}
				}else if(p.termEquals("Q", 0)){
					pair(data, active, affected, pairs);
				}else if(p.termEquals("T", 0)){
					require(p.terms()==6 && p.parseInt(4)==active.size() && p.parseInt(5)==pairs.size()
						&& data.contigs.add(((long)source<<32)|readID),
						"Every input/read must conserve all selected-feature and pair rows exactly once");
					for(int view=0; view<2; view++){
						data.counts[view][9]+=affected[view].size();
						affected[view].clear();
					}
					active.clear();
					pairs.clear();
				}else{throw new IllegalArgumentException("Unknown selected-geometry record type in "+path);}
			}
		}
		require(active.isEmpty(), "Selected geometry ended before its contig conservation row");
		return data;
	}

	/** Checks emitted pair arithmetic against retained P rows before counting it. */
	private void pair(final Data data, final ArrayList<Row> rows, final IntHashSet[] affected, final LongHashSet pairs){
		require(p.terms()==16, "Selected neighbour rows require the sixteen declared fields");
		final int ai=p.parseInt(4), bi=p.parseInt(5);
		require(ai>=0 && ai<rows.size() && bi>=0 && bi<rows.size() && ai!=bi, "Pair indices must identify distinct retained selected objects");
		final Row a=rows.get(ai), b=rows.get(bi);
		require(movedTarget(a) && pairs.add(((long)Math.min(ai, bi)<<32)|Math.max(ai, bi)),
			"Cobalamin overlap evidence requires a moved target and one row per unordered pair");
		final int edge=flag(6), same=flag(7), before=p.parseInt(9), after=p.parseInt(10), cap=p.parseInt(11);
		require(same==(a.strand==b.strand ? 1 : 0) && p.parseInt(8)==b.type && cap>=0
			&& before==overlap(a.start, a.stop, b.start, b.stop)
			&& after==overlap(a.reportStart, a.reportStop, b.reportStart, b.reportStop),
			"Saved overlap must describe these same selected objects in both coordinate states");
		final int oldOrder=p.parseInt(12), newOrder=p.parseInt(13), oldContains=p.parseInt(14), newContains=p.parseInt(15);
		require(oldOrder==order(a.start, a.stop, b.start, b.stop)
			&& newOrder==order(a.reportStart, a.reportStop, b.reportStart, b.reportStop)
			&& oldContains==containment(a.start, a.stop, b.start, b.stop)
			&& newContains==containment(a.reportStart, a.reportStop, b.reportStart, b.reportStop),
			"Saved order/containment must agree with the same retained P coordinates");
		for(int view=0; view<2; view++){
			if(view==1 && edge==0){continue;}
			final long[] c=data.counts[view];
			c[0]++;
			if(after>before){
				c[1]++;
				c[2]+=after-before;
				c[3]=Math.max(c[3], after-before);
				affected[view].add(ai);
				if(movedTarget(b)){affected[view].add(bi);}
			}
			if(before==0 && after>0){c[4]++;}
			if(before<=cap && after>cap){c[5]++;}
			if(oldOrder!=newOrder){c[6]++;}
			if(oldContains!=newContains){c[7]++;}
			if(same==0){c[8]++;}
		}
	}

	private boolean movedTarget(final Row row){
		final byte[] name=symbols.bytes(row.family);
		if(name.length<9){return false;}
		for(int i=0; i<9; i++){if(name[i]!=COBALAMIN[i]){return false;}}
		return (name.length==9 || name[9]=='_') && (row.start!=row.reportStart || row.stop!=row.reportStop);
	}
	private int flag(final int field){
		final int value=p.parseInt(field);
		require(value==0 || value==1, "Selection flags are literal0or1");
		return value;
	}
	private int symbol(final int field){
		p.setBounds(field);
		require(p.b()>p.a(), "Selected identity fields cannot be empty");
		return symbols.intern(p.line(), p.a(), p.b());
	}
	private static int overlap(final int a, final int b, final int c, final int d){return Math.max(0, Math.min(b, d)-Math.max(a, c)+1);}
	private static int order(final int a, final int b, final int c, final int d){return b==d ? Integer.compare(a, c) : Integer.compare(b, d);}
	private static int containment(final int a, final int b, final int c, final int d){return (a<=c && b>=d ? 1 : 0)|(c<=a && d>=b ? 2 : 0);}
	private static void metric(final Output output, final ByteBuilder b, final String name, final long value){output.writer.print(b.clear().append(name).tab().append(value).nl());}
	private static void writeStats(final Output output, final ByteBuilder b, final String prefix, final Data data){
		metric(output, b, prefix+"_moved_report_features", data.moved);
		for(int view=0; view<2; view++){
			for(int i=0; i<METRICS.length; i++){metric(output, b, prefix+(view==0 ? "_all_pairs_" : "_selected_edges_")+METRICS[i], data.counts[view][i]);}
		}
	}

	/** P rows include unchanged and unmatched calls; strand maps genomic ends to sense ends. */
	private void countMovement(final Data data, final Row row){
		final byte[] name=symbols.bytes(row.family);
		if(name.length<COBALAMIN.length){return;}
		for(int i=0; i<COBALAMIN.length; i++){if(name[i]!=COBALAMIN[i]){return;}}
		if(name.length>COBALAMIN.length && name[COBALAMIN.length]!='_'){return;}
		int member=-1;
		for(int i=0; i<MEMBER_BYTES.length; i++){if(Arrays.equals(name, MEMBER_BYTES[i])){member=i;break;}}
		require(member>=0, "Physical movement counts require an explicit registered cobalamin member, never inferred subtype");
		final boolean left=row.start!=row.reportStart;
		final boolean right=row.stop!=row.reportStop;
		final boolean five=row.strand==0 ? left : right;
		final boolean three=row.strand==0 ? right : left;
		for(int pass=0; pass<2; pass++){
			final long[] counts=data.movements[pass==0 ? member : 3];
			counts[0]++;
			if(five){counts[1]++;}
			if(three){counts[2]++;}
			if(five && three){counts[3]++;}
			if(five || three){counts[4]++;}
		}
	}

	private static void writeMovements(final Output output, final String view, final Data data){
		final ByteBuilder b=new ByteBuilder();
		for(int field=0; field<5; field++){
			require(data.movements[0][field]+data.movements[1][field]+data.movements[2][field]==data.movements[3][field],
				"Physical member counts must conserve the clan at every movement field");
		}
		for(int member=0; member<4; member++){
			final long[] counts=data.movements[member];
			require(counts[1]+counts[2]-counts[3]==counts[4] && counts[4]<=counts[0],
				"Single and both-end movements must conserve distinct selected calls");
			b.clear().append(view).tab().append(MEMBERS[member]);
			for(final long count:counts){b.tab().append(count);}
			output.writer.print(b.nl());
		}
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Data{
		final IntHashMap inputs=new IntHashMap();
		final LongHashSet contigs=new LongHashSet();
		final ArrayList<Row> rows=new ArrayList<Row>();
		final long[][] counts=new long[2][10];
		final long[][] movements=new long[4][5];
		long moved;
	}

	/** Retained output identity, not a transient wrapper for parsing a primitive. */
	private static final class Row implements Comparable<Row>{
		Row(final int source_, final long read_, final int contig_, final int strand_, final int type_, final int frame_,
				final int family_, final int model_, final int start_, final int stop_, final int reportStart_, final int reportStop_, final int score_){
			source=source_;
			read=read_;
			contig=contig_;
			strand=strand_;
			type=type_;
			frame=frame_;
			family=family_;
			model=model_;
			start=start_;
			stop=stop_;
			reportStart=reportStart_;
			reportStop=reportStop_;
			score=score_;
		}
		@Override
		public int compareTo(final Row other){
			int c=Integer.compare(source, other.source);
			if(c==0){c=Long.compare(read, other.read);}
			if(c==0){c=Integer.compare(contig, other.contig);}
			if(c==0){c=Integer.compare(strand, other.strand);}
			if(c==0){c=Integer.compare(type, other.type);}
			if(c==0){c=Integer.compare(frame, other.frame);}
			if(c==0){c=Integer.compare(family, other.family);}
			if(c==0){c=Integer.compare(model, other.model);}
			if(c==0){c=Integer.compare(start, other.start);}
			if(c==0){c=Integer.compare(stop, other.stop);}
			return c==0 ? Integer.compare(score, other.score) : c;
		}
		final int source, contig, strand, type, frame, family, model, start, stop, reportStart, reportStop, score;
		final long read;
	}

	private final Symbols symbols=new Symbols();
	private final boolean movements;
	private final LineParser1 p=new LineParser1('\t');
	private static final byte[] COBALAMIN="cobalamin".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	private static final String[] MEMBERS={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689", "cobalamin"};
	private static final byte[][] MEMBER_BYTES={MEMBERS[0].getBytes(java.nio.charset.StandardCharsets.US_ASCII),
		MEMBERS[1].getBytes(java.nio.charset.StandardCharsets.US_ASCII), MEMBERS[2].getBytes(java.nio.charset.StandardCharsets.US_ASCII)};
	private static final String[] METRICS={"pairs", "increased_pairs", "increased_bases", "max_increase_bases", "new_overlaps",
		"cap_crossings", "order_changes", "containment_changes", "opposite_strand_pairs", "affected_target_rnas"};
}
