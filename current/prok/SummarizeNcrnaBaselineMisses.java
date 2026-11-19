package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.IntObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Exclusive full-roster outcomes from completed baseline grading, without DNA
 * or verifier replay. Adapted from template.A_SampleBasic2 and the saved grader's
 * native byte prepass. Final matching takes precedence over pre-path membership;
 * then only a captured, overlapping window distinguishes the remaining classes.
 * Native truncation and incomplete coverage are independent subset views, never
 * additive populations. Rejected identity/HBM scores are not in these schemas
 * and are deliberately not invented from accepted-only or missing fields.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaBaselineMisses{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, SummarizeNcrnaBaselineMisses.class, false);
		try{
			final SummarizeNcrnaBaselineMisses tool=new SummarizeNcrnaBaselineMisses(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("NCRNA_BASELINE_MISS_CLASSES_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaBaselineMisses(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("population")){population=b;}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedtruncated")){expectedTruncated=Parse.parseIntKMG(b);}
			else if(a.equals("expectedincomplete")){expectedIncomplete=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown miss-report argument: "+arg);}
		}
		manifest=parser.in1;
		prefix=parser.out1;
		require(manifest!=null && prefix!=null && population!=null && expectedGenes>0
			&& expectedTruncated>=0 && expectedIncomplete>=0 && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=grade_manifest.tsv population=prepared_prefix out=fresh_prefix and explicit expectedgenes/truncated/incomplete without sampling");
		for(final String suffix:new String[]{".classes.tsv", ".genes.tsv.gz"}){
			require(!new File(prefix+suffix).exists(), "Miss reports require fresh output: "+prefix+suffix);
		}
		Shared.setThreads(1);
	}

	private void process(){
		loadOriginals();
		coordinateOrder=genes.toArray(new Gene[genes.size()]);
		Arrays.sort(coordinateOrder, (a, b)->{
			final int frame=compare(a.frame, b.frame);
			return frame!=0 ? frame : compare(a.contig, b.contig);
		});
		try(Input list=new Input(manifest); Output classes=new Output(prefix+".classes.tsv"); Output detail=new Output(prefix+".genes.tsv.gz")){
			RrnaResourceIO.header(list.reader, "library\tgrade_prefix\tcapture_manifest", manifest);
			classes.writer.println("library\tcohort\toriginal_genes\tno_captured_overlap\tcaptured_no_gene_acceptance\taccepted_not_final\tfinal_matched\tclass1_only_scheduled_rejections\tclass2_without_captured_containment\tfinal_without_unique_prepath\tclass2_qualifying_acceptance_not_unique\tclass2_accepted_overlap_fails_geometry\tclass2_other_or_unresolved_acceptance\tclass2_no_accepted_event\tmissed_final_overlap_fails_geometry\tmissed_final_qualifying_overlap");
			detail.writer.println("library\toriginal_gene\tnative_truncated\tincomplete_coverage\toutcome\tcaptured_overlap\tcaptured_containment\tscheduled_overlap\tunique_prepath_acceptance\tfinal_matched\taccepted_event_in_captured_overlap\taccepted_qualifying_geometry\taccepted_overlap_failing_geometry\tfinal_qualifying_geometry\tfinal_overlap_failing_geometry\treference_resolved");
			final ObjectSet<String> seen=new ObjectSet<String>(String.class);
			final LineParser1 entries=new LineParser1('\t');
			for(byte[] line=list.reader.nextLine(); line!=null; line=list.reader.nextLine()){
				entries.set(line);
				require(entries.terms()==3, "Each library requires its completed grade prefix and capture manifest");
				final String library=entries.parseString(0), grade=entries.parseString(1), captures=entries.parseString(2);
				require(library.matches("[A-Za-z0-9_]+") && seen.add(library) && new File(grade).isAbsolute(),
					"Library labels must be unique simple names and grade prefixes absolute");
				readGrade(grade);
				readEvents(captures);
				readAssociations(grade);
				readFinalCalls(grade);
				write(library, classes, detail);
			}
			require(!seen.isEmpty(), "At least one completed library grade is required");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Population       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadOriginals(){
		int truncated=0, incomplete=0;
		try(Input in=new Input(population+".original_calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, PrepareNcrnaRosterBaseline.ORIGINAL_HEADER+"\tbaseline_input_id\tbaseline_input_state", population);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==31 && p.termEquals("true", 7), "Use the complete accepted original-call population");
				final int id=p.parseInt(0), bits=p.parseInt(17);
				final boolean partial=bool(21);
				final Gene gene=new Gene(id, encoded(29), encoded(3), p.parseInt(4), p.parseInt(5), boolStrand(6), bits, partial);
				require(id>=0 && !ids.contains(id) && gene.start>0 && gene.stop>=gene.start && bits>=0,
					"Original genes need unique IDs and intact source coordinates");
				ids.put(id, gene);
				genes.add(gene);
				if(bits!=0){truncated++;}
				if(partial){incomplete++;}
			}
		}
		require(genes.size()==expectedGenes && truncated==expectedTruncated && incomplete==expectedIncomplete,
			"The full original, native-truncated and incomplete-coverage denominators must match their accepted counts");
	}

	private void readGrade(final String grade){
		for(final Gene gene:genes){gene.clear();}
		readSummary(grade);
		int count=0;
		try(Input in=new Input(grade+".genes.tsv.gz")){
			RrnaResourceIO.header(in.reader, SummarizeNcrnaRosterBaseline.GENES_HEADER, grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==26, "Use the complete original-gene baseline report");
				final Gene gene=ids.get(p.parseInt(0));
				require(gene!=null && !gene.seen && field(1, gene.frame) && field(2, gene.contig)
					&& p.parseInt(3)==gene.start && p.parseInt(4)==gene.stop && boolStrand(5)==gene.reverse
					&& p.parseInt(6)==gene.truncation, "A grade cannot omit, repeat or rebind an original gene");
				gene.seen=true;
				gene.accepted=p.termEquals("UNIQUE_BEST_ACCEPTANCE", 11);
				require(gene.accepted || p.termEquals("UNRESOLVED_INPUT", 11) || p.termEquals("UNRESOLVED_EXTRACTION", 11)
					|| p.termEquals("AMBIGUOUS_BEST_ACCEPTANCE", 11) || p.termEquals("NO_SCHEDULED_WINDOW", 11)
					|| p.termEquals("PARTIAL_WINDOW_NO_ACCEPTANCE", 11) || p.termEquals("WINDOW_WITHOUT_ACCEPTED_MODEL", 11),
					"Unknown proposal states cannot silently become verifier misses");
				gene.finalMatched=bool(21);
				gene.resolved=!p.termEquals("UNRESOLVED_INPUT", 11) && !p.termEquals("UNRESOLVED_EXTRACTION", 11);
				gene.scheduled=p.parseLong(9)>0;
				count++;
			}
		}
		require(count==expectedGenes, "Every original gene must appear once, even without a window or call");
	}

	private void readAssociations(final String grade){
		long count=0;
		try(Input in=new Input(grade+".gene_windows.tsv.gz")){
			RrnaResourceIO.header(in.reader, SummarizeNcrnaRosterBaseline.ASSOCIATION_HEADER, grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==9, "Use the complete native gene/window association report");
				final Gene gene=ids.get(p.parseInt(0));
				require(gene!=null && gene.seen && field(1, gene.frame), "Window associations must retain the same original gene and frame");
				final boolean captured=!p.termEquals("NA", 3), overlap=bool(4), contains=bool(5);
				require(!contains || overlap, "A containing positive-length window must overlap its original gene");
				if(captured && overlap){gene.captured=true;}
				if(captured && contains){gene.contained=true;}
				if(captured){
					p.setBounds(3);
					final Event event=findEvent(p.line(), p.a(), p.b());
					require(event!=null && Arrays.equals(event.frame, gene.frame) && event.reverse==gene.reverse,
						"Captured associations must bind their exact original event namespace, frame and strand");
					if(overlap && event.accepted){
						gene.acceptedEvent=true;
						if(gene.resolved && overlap(gene.start, gene.stop, event.start, event.stop)>0){
							if(qualifies(gene.start, gene.stop, event.start, event.stop)){gene.qualifyingEvent=true;}
							else{gene.geometryEvent=true;}
						}
					}
				}
				count++;
			}
		}
		require(count==expectedAssociations, "The association table must retain every row counted by the completed grade");
	}

	private void readSummary(final String grade){
		expectedAssociations=-1;
		expectedFinal=-1;
		expectedAccepted=-1;
		expectedCaptured=-1;
		expectedCalls=-1;
		boolean original=false;
		try(Input in=new Input(grade+".summary.tsv")){
			RrnaResourceIO.header(in.reader, "metric\tcount", grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==2, "Completed baseline summaries require named nonnegative counts");
				final long value=p.parseLong(1);
				require(value>=0, "A completed grade cannot contain negative counts");
				if(p.termEquals("original_genes", 0)){
					require(!original && value==expectedGenes, "The completed grade must retain the exact original population");
					original=true;
				}else if(p.termEquals("gene_window_associations_nonadditive_cost", 0)){
					require(expectedAssociations<0, "The completed grade cannot repeat its association count");
					expectedAssociations=value;
				}else if(p.termEquals("final_one_to_one_detected_genes", 0)){
					require(expectedFinal<0, "The completed grade cannot repeat its matched count");
					expectedFinal=value;
				}else if(p.termEquals("unique_best_members", 0)){
					require(expectedAccepted<0, "The completed grade cannot repeat its accepted count");
					expectedAccepted=value;
				}else if(p.termEquals("captured_window_events", 0)){
					require(expectedCaptured<0, "The completed grade cannot repeat its captured-event count");
					expectedCaptured=value;
				}else if(p.termEquals("final_call_rows", 0)){
					require(expectedCalls<0, "The completed grade cannot repeat its final-call count");
					expectedCalls=value;
				}
			}
		}
		require(original && expectedAssociations>=0 && expectedFinal>=0 && expectedAccepted>=0 && expectedCaptured>=0 && expectedCalls>=0,
			"Completed baseline counters are required to cross-check gene and association conservation");
	}

	/** Keep the captured namespace and real accepted OWN span. Sorting retained
	 * byte keys permits allocation-free joins from the association stream. */
	private void readEvents(final String manifest){
		final ArrayList<Event> order=new ArrayList<Event>();
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final LineParser1 entry=new LineParser1('\t');
		boolean unnamed=false;
		try(Input list=new Input(manifest)){
			RrnaResourceIO.header(list.reader, "capture_id\tcapture_prefix\treference_manifest", manifest);
			for(byte[] line=list.reader.nextLine(); line!=null; line=list.reader.nextLine()){
				entry.set(line);
				require(entry.terms()==3, "Use each library's exact saved capture manifest");
				final String name=entry.parseString(0), source=entry.parseString(1);
				require(names.add(name) && (name.isEmpty() || name.matches("[A-Za-z0-9_]+")) && new File(source).isAbsolute(),
					"Captured event namespaces must remain distinct and prefixes absolute");
				unnamed|=name.isEmpty();
				try(Input in=new Input(source+".results.tsv")){
					RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.RESULT_HEADER, source);
					for(byte[] event=in.reader.nextLine(); event!=null; event=in.reader.nextLine()){
						p.set(event);
						require(p.terms()==26, "Use the complete original completed-window result schema");
						if(p.termEquals("NA", 12)){continue;}
						field.clear();
						if(!name.isEmpty()){field.append(name).append(':');}
						p.setBounds(12);
						field.append(event, p.a(), p.b()-p.a());
						final byte[] id=field.toBytes(), frame=p.parseByteArray(1);
						final boolean accepted=p.termEquals("ACCEPTED_PROPOSAL", 13);
						final int strand=p.parseInt(5), start=accepted ? p.parseInt(22)+1 : -1, stop=accepted ? p.parseInt(23)+1 : -1;
						require((strand==0 || strand==1) && (!accepted || start>0 && stop>=start), "Accepted event geometry must retain explicit original coordinates");
						order.add(new Event(id, frame, strand==1, accepted, start, stop));
					}
				}
			}
		}
		require(!names.isEmpty() && (!unnamed || names.size()==1) && order.size()==expectedCaptured,
			"Every counted captured event must survive exactly once; an unnamed legacy capture must be alone");
		events=order.toArray(new Event[order.size()]);
		Arrays.sort(events, (a, b)->compare(a.id, b.id));
		for(int i=1; i<events.length; i++){require(compare(events[i-1].id, events[i].id)!=0, "Captured event identities cannot collide across inputs or shards");}
	}

	private Event findEvent(final byte[] row, final int from, final int to){
		int lo=0, hi=events.length;
		while(lo<hi){
			final int mid=(lo+hi)>>>1;
			final int cmp=compare(row, from, to, events[mid].id);
			if(cmp==0){return events[mid];}
			if(cmp<0){hi=mid;}else{lo=mid+1;}
		}
		return null;
	}

	/** An overlapping final call may fail reciprocal geometry or lose one-to-one
	 * matching. Neither is evidence that every model failed verification. */
	private void readFinalCalls(final String grade){
		long count=0;
		try(Input in=new Input(grade+".calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, SummarizeNcrnaRosterBaseline.CALLS_HEADER, grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==11, "Use the input-bound completed final-call grading table");
				p.setBounds(1);
				final int fa=p.a(), fb=p.b();
				p.setBounds(2);
				final int ca=p.a(), cb=p.b(), start=p.parseInt(3), stop=p.parseInt(4);
				final boolean reverse=boolStrand(5);
				int lo=0, hi=coordinateOrder.length;
				while(lo<hi){
					final int mid=(lo+hi)>>>1;
					if(compareTarget(line, fa, fb, ca, cb, coordinateOrder[mid])>0){lo=mid+1;}else{hi=mid;}
				}
				for(int i=lo; i<coordinateOrder.length && compareTarget(line, fa, fb, ca, cb, coordinateOrder[i])==0; i++){
					final Gene gene=coordinateOrder[i];
					if(!gene.resolved || gene.reverse!=reverse || overlap(gene.start, gene.stop, start, stop)==0){continue;}
					if(qualifies(gene.start, gene.stop, start, stop)){gene.finalQualifying=true;}
					else{gene.finalGeometry=true;}
				}
				count++;
			}
		}
		require(count==expectedCalls, "Final geometric diagnostics must retain every counted query call");
	}

	private static int compareTarget(byte[] row, int fa, int fb, int ca, int cb, Gene gene){
		final int frame=compare(row, fa, fb, gene.frame);
		return frame!=0 ? frame : compare(row, ca, cb, gene.contig);
	}
	private static int compare(final byte[] a, final byte[] b){return compare(a, 0, a.length, b);}
	private static int compare(final byte[] a, final int from, final int to, final byte[] b){
		final int length=Math.min(to-from, b.length);
		for(int i=0; i<length; i++){
			final int diff=(a[from+i]&255)-(b[i]&255);
			if(diff!=0){return diff;}
		}
		return Integer.compare(to-from, b.length);
	}
	private static int overlap(int a, int b, int c, int d){return Math.max(0, Math.min(b, d)-Math.max(a, c)+1);}
	private static boolean qualifies(int a, int b, int c, int d){
		final long shared=overlap(a, b, c, d);
		return 2*shared>=b-(long)a+1 && 2*shared>=d-(long)c+1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String library, final Output classes, final Output detail){
		final long[][] counts=new long[3][13];
		for(final Gene gene:genes){
			assert(gene.seen) : "All original genes must be joined before classifying misses";
			final int outcome=gene.finalMatched ? 3 : gene.accepted ? 2 : gene.captured ? 1 : 0;
			for(int cohort=0; cohort<3; cohort++){
				if(cohort==1 && gene.truncation==0 || cohort==2 && !gene.incomplete){continue;}
				final long[] cell=counts[cohort];
				cell[outcome]++;
				if(outcome==0 && gene.scheduled){cell[4]++;}
				if(outcome==1 && !gene.contained){cell[5]++;}
				if(gene.finalMatched && !gene.accepted){cell[6]++;}
				if(outcome==1){
					cell[gene.qualifyingEvent ? 7 : gene.geometryEvent ? 8 : gene.acceptedEvent ? 9 : 10]++;
				}
				if(!gene.finalMatched && gene.finalGeometry){cell[11]++;}
				if(!gene.finalMatched && gene.finalQualifying){cell[12]++;}
			}
			detail.writer.print(row.clear().append(library).tab().append(gene.id).tab().append(gene.truncation!=0).tab()
				.append(gene.incomplete).tab().append(OUTCOMES[outcome]).tab().append(gene.captured).tab()
				.append(gene.contained).tab().append(gene.scheduled).tab().append(gene.accepted).tab().append(gene.finalMatched).tab()
				.append(gene.acceptedEvent).tab().append(gene.qualifyingEvent).tab().append(gene.geometryEvent).tab()
				.append(gene.finalQualifying).tab().append(gene.finalGeometry).tab().append(gene.resolved).nl());
		}
		final int[] expected={expectedGenes, expectedTruncated, expectedIncomplete};
		require(counts[0][3]==expectedFinal && counts[0][2]+counts[0][3]-counts[0][6]==expectedAccepted,
			"Final matches and unique pre-path memberships must agree with their separately counted grade totals");
		for(int cohort=0; cohort<3; cohort++){
			final long[] cell=counts[cohort];
			require(cell[0]+cell[1]+cell[2]+cell[3]==expected[cohort], "The four exclusive classes must conserve each declared cohort independently");
			require(cell[7]+cell[8]+cell[9]+cell[10]==cell[1], "The separate class2 geometry/other-acceptance/no-acceptance subclasses must conserve class2");
			row.clear().append(library).tab().append(COHORTS[cohort]).tab().append(expected[cohort]);
			for(final long count:cell){row.tab().append(count);}
			classes.writer.print(row.nl());
		}
	}

	private boolean bool(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Saved boolean fields must be explicit true/false");
		return p.termEquals("true", column);
	}
	private boolean boolStrand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Original strands cannot be inferred from absent data");
		return p.termEquals("-", column);
	}
	private byte[] encoded(final int column){
		p.setBounds(column);
		return BuildNcrnaReference.field(field.clear(), p.line(), p.a(), p.b()).toBytes();
	}
	private boolean field(final int column, final byte[] expected){
		p.setBounds(column);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Gene{
		Gene(int id_, byte[] frame_, byte[] contig_, int start_, int stop_, boolean reverse_, int truncation_, boolean incomplete_){
			id=id_;
			frame=frame_;
			contig=contig_;
			start=start_;
			stop=stop_;
			reverse=reverse_;
			truncation=truncation_;
			incomplete=incomplete_;
		}
		void clear(){
			seen=false;
			accepted=false;
			finalMatched=false;
			captured=false;
			contained=false;
			scheduled=false;
			resolved=false;
			acceptedEvent=false;
			qualifyingEvent=false;
			geometryEvent=false;
			finalQualifying=false;
			finalGeometry=false;
		}
		final int id, start, stop, truncation;
		final byte[] frame, contig;
		final boolean reverse, incomplete;
		boolean seen, accepted, finalMatched, captured, contained, scheduled;
		boolean resolved, acceptedEvent, qualifyingEvent, geometryEvent, finalQualifying, finalGeometry;
	}
	private static final class Event{
		Event(byte[] id_, byte[] frame_, boolean reverse_, boolean accepted_, int start_, int stop_){
			id=id_;
			frame=frame_;
			reverse=reverse_;
			accepted=accepted_;
			start=start_;
			stop=stop_;
		}
		final byte[] id, frame;
		final boolean reverse, accepted;
		final int start, stop;
	}
	private String manifest, population, prefix;
	private int expectedGenes, expectedTruncated=-1, expectedIncomplete=-1;
	private long expectedAssociations, expectedFinal, expectedAccepted;
	private long expectedCaptured, expectedCalls;
	private Event[] events;
	private Gene[] coordinateOrder;
	private final IntObjectMap<Gene> ids=new IntObjectMap<Gene>(16, Gene.class);
	private final ArrayList<Gene> genes=new ArrayList<Gene>();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder(), field=new ByteBuilder();
	private static final String[] COHORTS={"ALL_ORIGINALS", "NATIVE_TRUNCATED", "INCOMPLETE_COVERAGE"};
	private static final String[] OUTCOMES={"NO_CAPTURED_OVERLAP", "CAPTURED_NO_GENE_ACCEPTANCE", "ACCEPTED_NOT_FINAL", "FINAL_MATCHED"};
}
