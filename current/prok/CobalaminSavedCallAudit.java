package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;

import map.IntObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Audits saved normalized cobalamin calls without calling or grading again.
 * Forked from template.A_SampleByteFile; one metadata thread owns native readers,
 * retained call records, primitive counters and output builders. Exact identity
 * is source/contig/strand/coordinates, independently of member/model labels.
 * Spatial links mean at least one shared base on that same source and strand.
 * They describe final-call correspondence, never prove candidate ancestry.
 * Only unique bidirectional spatial counterparts get masked/matched provenance;
 * mixed, multiple, prior-unmatched and unchanged-coordinate grade changes remain
 * explicitly separate. No existing match, mask or reference assignment changes.
 * Legacy unknown scored eligibility stays NA; classifications use saved grades
 * and physical coordinates, never inferred eligibility.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSavedCallAudit{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, CobalaminSavedCallAudit.class, false);
		try{
			final CobalaminSavedCallAudit tool=new CobalaminSavedCallAudit(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("SAVED_CALL_AUDIT_PASS before="+tool.calls[0].length+" after="+tool.calls[1].length+" time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CobalaminSavedCallAudit(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("expectedinputs")){
				expectedInputs=Parse.parseIntKMG(b);
			}else{require(parser.parse(arg, a, b), "Unknown saved-call audit option: "+arg);}
		}
		prefixes[0]=parser.in1;
		prefixes[1]=parser.in2;
		directory=parser.out1;
		require(prefixes[0]!=null && prefixes[1]!=null && directory!=null && expectedInputs>0
			&& !parser.overwrite && !parser.append && parser.out2==null && parser.maxReads<0 && parser.samplerate==1
			&& !new File(directory).exists(),
			"Supply complete in1=OFF_prefix in2=ON_prefix expectedinputs=N out=fresh_directory without filtering or overwrite");
		Shared.setThreads(1);
	}

	private void process(){
		readSources(0);
		readSources(1);
		for(int side=0; side<2; side++){
			calls[side]=readCalls(side);
			Arrays.sort(calls[side], ORDER);
			for(int i=1; i<calls[side].length; i++){
				require(!exact(calls[side][i-1], calls[side][i]), "Duplicate coordinates in normalized side "+side);
			}
		}
		for(final Source source:sources){
			for(int side=0; side<2; side++){
				require(Arrays.equals(source.expected[side], source.actual[side]),
					"Saved query outcomes do not conserve declared genome counts: "+new String(symbols.bytes(source.id)));
			}
		}
		require(new File(directory).mkdir(), "Cannot create fresh audit directory: "+directory);
		try(Output pairs=output("overlaps.tsv"); Output on=output("on_calls.tsv");
				Output unmatched=output("on_unmatched.tsv"); Output off=output("off_calls.tsv")){
			pairs.writer.println("on_row\toff_row\tinput_id\tseqid\tstrand\ton_start\ton_stop\toff_start\toff_stop\ton_member\toff_member\ton_outcome\toff_outcome\tintersection_bp\tfraction_on\tfraction_off\texact_coordinates");
			link(pairs);
			final String header="query_row\tinput_id\tseqid\tstrand\tstart\tstop\tmember\tmodel\toutcome\tscored_eligible\tmatched_reference_id\texact_other_row\texact_other_member\texact_other_outcome\toverlap_matched\toverlap_masked\toverlap_unmatched\tclassification";
			on.writer.println(header);
			unmatched.writer.println(header);
			off.writer.println(header);
			for(final Call call:calls[1]){
				call.category=category(call);
				categories[call.member][call.category]++;
				callRow(call, true, on);
				if(call.grade()==2){callRow(call, true, unmatched);}
			}
			for(final Call call:calls[0]){callRow(call, false, off);}
		}
		writeCounts();
		writeSpatialCounts();
		try(Output complete=output("complete.tsv")){
			complete.writer.println("state\tbefore_calls\tafter_calls\tinputs\tidentity\toverlap");
			complete.writer.print(buffer.clear().append("ALL_CALLS_AND_UNMATCHED_PARTITIONS_CONSERVED\t").append(calls[0].length)
				.tab().append(calls[1].length).tab().append(sources.size())
				.append("\tSOURCE_CONTIG_STRAND_COORDINATES_IGNORE_MEMBER_MODEL\tSAME_SOURCE_CONTIG_STRAND_AT_LEAST_ONE_BASE_NOT_CAUSAL\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Complete Inputs       ----------------*/
	/*--------------------------------------------------------------*/

	/** Both normalized runs must contain the same complete input roster, including
	 * inputs with no calls. Their independently declared grade counts are checked.
	 */
	private void readSources(final int side){
		int count=0;
		try(Input input=new Input(prefixes[side]+".genomes.tsv")){
			RrnaResourceIO.header(input.reader, GENOMES, prefixes[side]);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				parse(row, 5);
				final int id=symbol(0);
				Source source=bySource.get(id);
				if(side==0){
					require(source==null, "Repeated OFF input identity");
					source=new Source(id);
					bySource.put(id, source);
					sources.add(source);
				}else{require(source!=null && !source.afterSeen, "ON input is absent from OFF or duplicated");}
				for(int i=0; i<4; i++){source.expected[side][i]=number(i+1);}
				require(source.expected[side][0]==source.expected[side][1]+source.expected[side][2]+source.expected[side][3],
					"Declared calls must partition into matched, masked and unmatched");
				if(side==1){source.afterSeen=true;}
				count++;
			}
		}
		require(count==expectedInputs, "Saved input roster count differs from expectedinputs: "+count);
	}

	/** Retains one small record per final call, with immutable interned metadata.
	 * Row numbers identify source TSV lines, not biological IDs or cross-run identity.
	 */
	private Call[] readCalls(final int side){
		final ArrayList<Call> result=new ArrayList<Call>();
		try(Input input=new Input(prefixes[side]+".queries.tsv")){
			RrnaResourceIO.header(input.reader, QUERIES, prefixes[side]);
			int line=1;
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				line++;
				parse(row, 10);
				final Source source=bySource.get(symbol(0));
				require(source!=null, "Saved call names an undeclared input at query row "+line);
				final int seqid=symbol(1), start=number(2), stop=number(3);
				final char strand=p.termEquals("+", 4) ? '+' : p.termEquals("-", 4) ? '-' : '.';
				final int member=member(5), model=symbol(6), outcome=outcome(7);
				require(start>0 && stop>=start && strand!='.', "Invalid saved call coordinates/strand at query row "+line);
				require(p.termEquals("true", 8) || p.termEquals("false", 8) || p.termEquals("NA", 8),
					"Invalid scored eligibility at query row "+line+"; require true, false or explicit NA");
				final int eligible=p.termEquals("NA", 8) ? -1 : p.termEquals("true", 8) ? 1 : 0;
				final boolean named=!(p.termEquals("NA", 9) || p.termEquals(".", 9));
				require((outcome==0)==named, "Only an already graded TP may carry a reference assignment");
				final Call call=new Call(line, source, seqid, start, stop, strand, member, model, outcome,
					eligible, symbol(9));
				result.add(call);
				source.actual[side][0]++;
				source.actual[side][1+call.grade()]++;
				totals[side][member][0]++;
				totals[side][member][1+call.grade()]++;
			}
		}
		return result.toArray(new Call[result.size()]);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Spatial Correspondence   ----------------*/
	/*--------------------------------------------------------------*/

	/** Joins source/contig/strand groups sorted by start. Each group is small;
	 * only baseline intervals starting before the ON stop can overlap. All links
	 * are retained, including one-to-many and mixed-grade links. No truth is rematched.
	 */
	private void link(final Output output){
		final Call[] before=calls[0];
		int from=0, to=0;
		Call prior=null;
		for(final Call on:calls[1]){
			if(prior==null || group(prior, on)!=0){
				while(from<before.length && group(before[from], on)<0){from++;}
				to=from;
				while(to<before.length && group(before[to], on)==0){to++;}
				prior=on;
			}
			for(int i=from; i<to && before[i].start<=on.stop; i++){
				final Call off=before[i];
				final int bp=Math.min(off.stop, on.stop)-Math.max(off.start, on.start)+1;
				if(bp<=0){continue;}
				on.overlaps[off.grade()]++;
				off.overlaps[on.grade()]++;
				if(off.overlapCount()>1){
					assert(off.only!=null) : "A second ON overlap requires the retained first counterpart";
					on.sharesOff=true;
					off.only.sharesOff=true;
				}
				on.only=off;
				off.only=on;
				if(exact(off, on)){
					require(off.exact==null && on.exact==null, "Exact coordinates require unique final calls on each side");
					off.exact=on;
					on.exact=off;
				}
				buffer.clear().append(on.row).tab().append(off.row).tab().append(symbols.bytes(on.source.id)).tab()
					.append(symbols.bytes(on.seqid)).tab().append(on.strand).tab().append(on.start).tab().append(on.stop)
					.tab().append(off.start).tab().append(off.stop).tab().append(MEMBERS[on.member]).tab().append(MEMBERS[off.member])
					.tab().append(OUTCOMES[on.outcome]).tab().append(OUTCOMES[off.outcome]).tab().append(bp)
					.tab().append(bp/(double)(on.stop-on.start+1), 8).tab().append(bp/(double)(off.stop-off.start+1), 8)
					.tab().append(exact(off, on) ? "true" : "false").nl();
				output.writer.print(buffer);
			}
		}
	}

	/** Disjoint ON-unmatched categories. Exact retained unmatched takes precedence;
	 * an unchanged span with a different grade is not called a moved endpoint.
	 */
	private static int category(final Call on){
		assert(on.source!=null) : "Every audit classification must retain its validated original input";
		if(on.grade()!=2){return 0;}
		if(on.exact!=null){return on.exact.grade()==2 ? 1 : 7;}
		if(on.overlapCount()==0){return 2;}
		if(on.overlapCount()!=1 || on.only.overlapCount()!=1){return 6;}
		if(on.only.grade()==1){return 3;}
		if(on.only.grade()==0){return 4;}
		return 5;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	private void callRow(final Call call, final boolean after, final Output output){
		assert(call.exact==null || exact(call, call.exact)) : "An exact counterpart must preserve source/contig/strand/coordinates despite member/model changes";
		buffer.clear().append(call.row).tab().append(symbols.bytes(call.source.id)).tab().append(symbols.bytes(call.seqid))
			.tab().append(call.strand).tab().append(call.start).tab().append(call.stop).tab().append(MEMBERS[call.member])
			.tab().append(symbols.bytes(call.model)).tab().append(OUTCOMES[call.outcome])
			.tab().append(call.eligible<0 ? "NA" : call.eligible==1 ? "true" : "false").tab().append(symbols.bytes(call.reference)).tab();
		if(call.exact==null){buffer.append("NA\tNA\tNA");}
		else{buffer.append(call.exact.row).tab().append(MEMBERS[call.exact.member]).tab().append(OUTCOMES[call.exact.outcome]);}
		for(final int count:call.overlaps){buffer.tab().append(count);}
		buffer.tab().append(after ? CATEGORIES[call.category] : offState(call)).nl();
		output.writer.print(buffer);
	}

	/** OFF-only unmatched is defined by the same exact identity as ON additions.
	 * Changed-coordinate overlaps stay visible separately from actual absence.
	 */
	private static String offState(final Call off){
		if(off.grade()!=2){return "NOT_UNMATCHED";}
		if(off.exact!=null){return off.exact.grade()==2 ? "RETAINED_EXACT_UNMATCHED" : "EXACT_OTHER_GRADE";}
		return off.overlapCount()==0 ? "NO_ON_CALL_AT_LOCUS" : "CHANGED_COORDINATE_ON_OVERLAP";
	}

	/** Conservation is checked independently for each emitted member and globally.
	 * Retained exact calls can change labels; incoming/outgoing switches explicitly
	 * reconcile the per-member net instead of making such calls look newly born.
	 */
	private void writeCounts(){
		final long[] beforeRetained=new long[3], afterRetained=new long[3];
		final long[] labelIn=new long[3], labelOut=new long[3];
		for(final Call off:calls[0]){
			if(off.grade()==2 && off.exact!=null && off.exact.grade()==2){
				beforeRetained[off.member]++;
				afterRetained[off.exact.member]++;
				if(off.member!=off.exact.member){
					labelOut[off.member]++;
					labelIn[off.exact.member]++;
				}
			}
		}
		try(Output members=output("members.tsv"); Output bins=output("categories.tsv")){
			members.writer.println("member\toff_calls\ton_calls\toff_matched\ton_matched\toff_masked\ton_masked\toff_unmatched\ton_unmatched\tretained_unmatched_off\tretained_unmatched_on\ton_only_unmatched\toff_only_unmatched\tlabel_switch_in\tlabel_switch_out\tunmatched_net_change");
			bins.writer.println("member\tclassification\tcalls");
			final long[] all=new long[15], allCategories=new long[CATEGORIES.length];
			for(int member=0; member<3; member++){
				final long[] a=totals[0][member], z=totals[1][member], c=categories[member];
				long unmatched=0, added=0;
				for(int i=1; i<c.length; i++){
					unmatched+=c[i];
					if(i!=1){added+=c[i];}
				}
				final long removed=a[3]-beforeRetained[member];
				require(a[0]==a[1]+a[2]+a[3] && z[0]==z[1]+z[2]+z[3]
					&& unmatched==z[3] && c[0]==z[1]+z[2] && c[1]==afterRetained[member]
					&& z[3]-a[3]==added-removed+labelIn[member]-labelOut[member],
					"Call categories or exact-identity additions/departures do not conserve member "+MEMBERS[member]);
				final long[] values={a[0], z[0], a[1], z[1], a[2], z[2], a[3], z[3],
					beforeRetained[member], afterRetained[member], added, removed,
					labelIn[member], labelOut[member], z[3]-a[3]};
				countRow(members, MEMBERS[member], values);
				for(int i=0; i<all.length; i++){all[i]+=values[i];}
				for(int i=0; i<c.length; i++){
					bins.writer.print(buffer.clear().append(MEMBERS[member]).tab().append(CATEGORIES[i]).tab().append(c[i]).nl());
					allCategories[i]+=c[i];
				}
			}
			require(all[0]==calls[0].length && all[1]==calls[1].length && all[8]==all[9]
				&& all[14]==all[10]-all[11], "Global unmatched additions/departures and label switches must reconcile exactly");
			countRow(members, "ALL", all);
			for(int i=0; i<allCategories.length; i++){
				bins.writer.print(buffer.clear().append("ALL").tab().append(CATEGORIES[i]).tab().append(allCategories[i]).nl());
			}
		}
		outstream.println("Final-call overlap is descriptive, not causal ancestry; exact identity ignores member/model labels.");
	}

	private void countRow(final Output output, final String label, final long[] values){
		buffer.clear().append(label);
		for(final long value:values){buffer.tab().append(value);}
		output.writer.print(buffer.nl());
	}

	/** Counts geometric sharing, not ancestry: at least two ON calls overlap one OFF call. */
	private void writeSpatialCounts(){
		final long[][] counts=new long[3][6];
		for(final Call on:calls[1]){
			if(on.grade()!=2){continue;}
			final long[] row=counts[on.member];
			row[0]++;
			if(on.exact!=null && on.exact.grade()==2){row[1]++;}
			else{
				row[2]++;
				row[on.overlapCount()==0 ? 3 : 4]++;
				if(on.sharesOff){row[5]++;}
			}
		}
		try(Output output=output("unmatched_spatial.tsv")){
			output.writer.println("member\ton_unmatched\tretained_exact_unmatched\ton_only_unmatched\ton_only_no_off_overlap\ton_only_off_overlap\ton_only_shares_off_with_another_on");
			final long[] total=new long[6];
			for(int member=0; member<counts.length; member++){
				final long[] row=counts[member];
				require(row[0]==row[1]+row[2] && row[2]==row[3]+row[4] && row[5]<=row[4]
					&& row[0]==totals[1][member][3], "Spatial headline must conserve every saved unmatched call for "+MEMBERS[member]);
				countRow(output, MEMBERS[member], row);
				for(int i=0; i<row.length; i++){total[i]+=row[i];}
			}
			countRow(output, "ALL", total);
		}
	}
	private Output output(final String name){return new Output(new File(directory, name).getPath());}
	private void parse(final byte[] row, final int columns){
		p.set(row);
		require(p.terms()==columns, "Normalized table schema changed");
	}
	private int number(final int term){
		p.setBounds(term);
		return NcrnaNumbers.unsignedInt(p.line(), p.a(), p.b());
	}
	private int symbol(final int term){
		p.setBounds(term);
		require(p.b()>p.a(), "Empty saved identity field");
		return symbols.intern(p.line(), p.a(), p.b());
	}
	private int member(final int term){
		for(int i=0; i<MEMBERS.length; i++){if(p.termEquals(MEMBERS[i], term)){return i;}}
		throw new IllegalArgumentException("Unknown emitted cobalamin member");
	}
	private int outcome(final int term){
		for(int i=0; i<OUTCOMES.length; i++){if(p.termEquals(OUTCOMES[i], term)){return i;}}
		throw new IllegalArgumentException("Unknown normalized grade outcome");
	}
	private static int group(final Call a, final Call b){
		int result=Integer.compare(a.source.id, b.source.id);
		if(result==0){result=Integer.compare(a.seqid, b.seqid);}
		return result==0 ? Character.compare(a.strand, b.strand) : result;
	}
	private static boolean exact(final Call a, final Call b){return group(a, b)==0 && a.start==b.start && a.stop==b.stop;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Retained Records      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final int id_){id=id_;}
		final int id;
		final long[][] expected=new long[2][4], actual=new long[2][4];
		boolean afterSeen;
	}
	private static final class Call{
		Call(final int row_, final Source source_, final int seqid_, final int start_, final int stop_,
				final char strand_, final int member_, final int model_, final int outcome_, final int eligible_, final int reference_){
			row=row_;
			source=source_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			member=member_;
			model=model_;
			outcome=outcome_;
			eligible=eligible_;
			reference=reference_;
		}
		int grade(){return Math.min(2, outcome);}
		int overlapCount(){return overlaps[0]+overlaps[1]+overlaps[2];}
		final Source source;
		final int row, seqid, start, stop, member, model, outcome, reference;
		final char strand;
		/** Original eligibility: -1 unknown, 0 false, 1 true. It never selects a call. */
		final int eligible;
		/** Matched/masked/unmatched counts of opposite-side physical overlaps. */
		final int[] overlaps=new int[3];
		Call exact, only;
		/** An OFF counterpart also overlaps another ON call; no causal claim. */
		boolean sharesOff;
		int category;
	}
	private final String[] prefixes=new String[2];
	private String directory;
	private int expectedInputs=-1;
	private final Symbols symbols=new Symbols();
	private final IntObjectMap<Source> bySource=new IntObjectMap<Source>(Source.class);
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final Call[][] calls=new Call[2][];
	/** Side/member/total-matched-masked-unmatched. */
	private final long[][][] totals=new long[2][3][4];
	private final long[][] categories=new long[3][8];
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final PrintStream outstream;
	private static final Comparator<Call> ORDER=(a, b)->{
		int result=group(a, b);
		if(result==0){result=Integer.compare(a.start, b.start);}
		return result==0 ? Integer.compare(a.stop, b.stop) : result;
	};
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] OUTCOMES={"TP", "MASKED", "UNMATCHED_OTHER", "UNMATCHED_DUPLICATE"};
	private static final String[] CATEGORIES={"NOT_UNMATCHED", "RETAINED_EXACT_UNMATCHED", "A_NO_OFF_CALL_AT_LOCUS",
		"B_UNIQUE_OFF_MASKED_OVERLAP", "C_UNIQUE_OFF_MATCHED_OVERLAP", "D_PRIOR_UNMATCHED_SHIFTED",
		"D_MULTIPLE_OR_MIXED_COUNTERPARTS", "D_EXACT_COORDINATES_CHANGED_GRADE"};
	private static final String GENOMES="input_id\tcalls\tmatched\tmasked\tunmatched";
	private static final String QUERIES="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id";
}
