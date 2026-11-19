package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ReadWrite;
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

/** Makes nested panels from an existing all-genome permutation, retaining negative
 * inputs and source-qualified raw producing-CM truth, including native truncations.
 * Also projects every admitted original ID into the legacy seed-coverage ledger.
 * Reads metadata only: no genome sequence, extraction, annotation or truth adjudication.
 * Lifecycle and argument parsing follow template.A_SampleByteFile; native text owners
 * close and check every reader/writer. This tool has one metadata consumer thread.
 * @author Raiden
 */
public final class CobalaminNativePanels{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, CobalaminNativePanels.class, false);
		try{
			final CobalaminNativePanels tool=new CobalaminNativePanels(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("COBALAMIN_NATIVE_PANELS_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private CobalaminNativePanels(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("manifest")){
				manifest=b;
			}else if(a.equals("counts")){
				countsPath=b;
			}else if(a.equals("loci")){
				lociPath=b;
			}else if(a.equals("sourceids")){
				sourcePath=b;
			}else if(a.equals("expectedgenomes")){
				expectedGenomes=Parse.parseIntKMG(b);
			}else if(a.equals("expected")){
				expected=parseCounts(b, 3);
			}else if(a.equals("rawexpected")){
				rawExpected=parseCounts(b, 3);
			}else if(a.equals("panels")){
				panelSizes=parseCounts(b, -1);
			}else if(a.equals("chunkstart")){
				chunkStart=Parse.parseIntKMG(b);
			}else if(a.equals("chunksize")){
				chunkSize=Parse.parseIntKMG(b);
			}else if(a.equals("completedprefix")){
				completedPrefix=b;
			}else{require(parser.parse(arg, a, b), "Unknown panel option: "+arg);}
		}
		orderPath=parser.in1;
		directory=parser.out1;
		require(orderPath!=null && manifest!=null && countsPath!=null && lociPath!=null && sourcePath!=null
			&& directory!=null && expectedGenomes>0 && expected!=null && rawExpected!=null && panelSizes!=null,
			"Bind in=saved_permutation manifest= counts= loci= sourceids= out=fresh_directory expectedgenomes= expected=A,B,C rawexpected=A,B,C panels=N,M");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Metadata populations cannot be overwritten, appended or sampled again");
		int previous=0;
		for(final int size:panelSizes){
			require(size>previous && size<=expectedGenomes, "Nested panel sizes must increase within the complete roster");
			previous=size;
		}
		for(int m=0; m<3; m++){require(expected[m]>=0 && rawExpected[m]>=expected[m], "Raw counts include every admitted original");}
		require((chunkSize==0 && chunkStart==0 && completedPrefix==null)
			|| (chunkSize>0 && chunkStart>0 && chunkStart<expectedGenomes && completedPrefix!=null),
			"Optional suffix chunks require chunksize>0, 0<chunkstart<expectedgenomes and completedprefix=accepted_sources");
		require(!new File(directory).exists(), "Outputs must be fresh: "+directory);
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_PIGZ=false;
		ReadWrite.setZipThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		loadOrder();
		loadManifest();
		loadCounts();
		if(chunkSize>0){validateCompletedPrefix();}
		require(new File(directory).mkdir(), "Could not create fresh metadata directory: "+directory);
		writeSources();
		writeTruth();
		writePanels();
		if(chunkSize>0){writeChunks();}
		try(Output policy=output("POLICY.txt")){
			policy.writer.println("Panels are prefixes of the saved full-roster permutation; no new sampling or positivity/taxonomy filter.");
			policy.writer.println("Protection truth is operational raw producing-CM member coordinates, including native truncations, not independent biological subtype truth.");
			policy.writer.println("Scored/dont-care grading remains the separate canonical NEW-GFF union. No member identity is inferred from reference_model; per-member scored-positive counts are unavailable.");
			policy.writer.println("Seed ledger TRAIN means every original is a resource-construction input. Run NcrnaSeedCoverage fold=ALL representatives=f; there is no held-out subset or sequence deduplication. Strata are producing-CM members.");
		}
		// Completion is last: partial directories never masquerade as accepted panels.
		try(Output complete=output("complete.tsv")){
			complete.writer.println("genomes\tadmitted_originals\traw_member_loci\tstate");
			complete.writer.print(buffer.clear().append(assemblies.size()).tab().append(sum(expected)).tab()
				.append(sum(rawExpected)).append("\tMETADATA_ONLY_GRADING_PENDING\n"));
		}
		outstream.println("Prepared nested metadata panels for "+assemblies.size()+" genomes, retaining zero-locus inputs.");
	}

	/** Retains the supplied order exactly; this method never shuffles or opens FASTAs. */
	private void loadOrder(){
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		try(Input in=new Input(orderPath)){
			RrnaResourceIO.header(in.reader, ORDER_HEADER, orderPath);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==6, "Saved all-genome order requires six fields");
				final int input=p.parseInt(1), admitted=p.parseInt(3);
				final String id=p.parseString(0);
				require(input>=0 && admitted>=0 && !byInput.contains(input) && names.add(id), "Repeated or invalid assembly/input in saved permutation: "+id);
				final Assembly assembly=new Assembly(id, input, assemblies.size(), admitted, p.parseString(4), p.parseString(5));
				require(assembly.fasta.startsWith("/") && assembly.source.startsWith("/"), "Panel paths must retain absolute source bindings");
				assemblies.add(assembly);
				byInput.put(input, assembly);
			}
		}
		require(assemblies.size()==expectedGenomes, "The saved order must contain every selected genome: "+assemblies.size());
	}

	/** Binds canonical reference paths without interpreting their dispositions or subtypes. */
	private void loadManifest(){
		int count=0;
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, manifest);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==12, "Selected manifest requires twelve fields");
				final Assembly assembly=assembly(1, 0);
				require(assembly.gff==null && p.termEquals(assembly.fasta, 3) && p.termEquals(assembly.source, 7), "Selected manifest and saved permutation disagree: "+assembly.id);
				assembly.gff=p.parseString(4);
				require(assembly.gff.startsWith("/"), "Canonical NEW-GFF path must be absolute: "+assembly.id);
				count++;
			}
		}
		require(count==expectedGenomes, "Reference bindings must cover the complete selected genome roster");
	}

	/** Uses saved zero-inclusive counts, then checks them against both locus ledgers. */
	private void loadCounts(){
		final long[] raw=new long[3], eligible=new long[3];
		try(Input in=new Input(countsPath)){
			RrnaResourceIO.header(in.reader, COUNTS_HEADER, countsPath);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==8, "Saved per-genome member counts require eight fields");
				final Assembly assembly=assembly(1, 0);
				final int m=member(3), nativeCount=p.parseInt(5), admitted=p.parseInt(6), truncated=p.parseInt(7);
				require((assembly.mask&(1<<m))==0 && nativeCount>=0 && admitted>=0 && truncated>=0
					&& (long)admitted+truncated==nativeCount && p.parseLong(4)==nativeCount,
					"Each genome/member needs one conserved native/eligible/truncation count: "+assembly.id);
				assembly.mask|=1<<m;
				assembly.raw[m]=nativeCount;
				assembly.eligible[m]=admitted;
				raw[m]+=nativeCount;
				eligible[m]+=admitted;
			}
		}
		for(final Assembly assembly:assemblies){
			require(assembly.mask==7 && sum(assembly.eligible)==assembly.admitted, "Zero-inclusive member counts must cover every saved input: "+assembly.id);
		}
		for(int m=0; m<3; m++){
			require(raw[m]==rawExpected[m] && eligible[m]==expected[m], "Declared member totals changed for "+MEMBERS[m]);
		}
	}

	/** Keeps every original ID as its own representative; equal DNA is never examined. */
	private void writeSources(){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(sourcePath); Output ledger=output("seed_coverage.tsv")){
			RrnaResourceIO.header(in.reader, SOURCE_HEADER, sourcePath);
			ledger.writer.println("locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==11, "Original source-ID provenance requires eleven fields");
				final Assembly assembly=assembly(3, 2);
				final int m=member(1), start=p.parseInt(5), stop=p.parseInt(6);
				final String id=p.parseString(0);
				require(id.startsWith(MEMBERS[m]+"_input_"+assembly.input+"_locus_") && ids.add(id)
					&& start>0 && stop>=start && (long)stop-start+1==p.parseInt(8)
					&& (p.termEquals("+", 7) || p.termEquals("-", 7)), "Original ID/geometry must remain unique and source-bound: "+id);
				assembly.originals[m]++;
				ledger.writer.print(buffer.clear().append(id).tab().append(assembly.id).tab().append(assembly.id).tab()
					.append(MEMBERS[m]).tab().append(id).append("\tTRAIN\tfalse\n"));
			}
		}
		for(final Assembly assembly:assemblies){
			for(int m=0; m<3; m++){
				require(assembly.originals[m]==assembly.eligible[m], "Source-ID ledger lost an admitted original: "+assembly.id+" "+MEMBERS[m]);
			}
		}
	}

	/** Streams the complete native ledger once into all nested truth views. */
	private void writeTruth(){
		try(NcrnaTextIO.Outputs outputs=new NcrnaTextIO.Outputs(panelSizes.length); Input in=new Input(lociPath)){
			final fileIO.ByteStreamWriter[] writers=new fileIO.ByteStreamWriter[panelSizes.length];
			for(int i=0; i<writers.length; i++){
				writers[i]=outputs.open(fileIO.FileFormat.testOutput(path("panel"+panelSizes[i]+".raw_native.truth.tsv"), fileIO.FileFormat.TXT, null, true, false, false, false));
				writers[i].println("seqid\tfamily\tstart\tstop\tstrand\tinput_ordinal0");
			}
			RrnaResourceIO.header(in.reader, LOCUS_HEADER, lociPath);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==13, "Raw native coordinates require thirteen fields");
				final Assembly assembly=assembly(1, 0);
				final int m=member(4), start=p.parseInt(6), stop=p.parseInt(7), full=p.parseInt(10);
				require(p.termEquals(assembly.source, 3) && start>0 && stop>=start && p.parseInt(9)==1
					&& (p.termEquals("+", 8) || p.termEquals("-", 8))
					&& (p.termEquals("true", 11) ? full==1 : p.termEquals("false", 11) && full==0), "Raw native source/coordinates/admission changed: "+assembly.id);
				assembly.observedRaw[m]++;
				assembly.observedEligible[m]+=full;
				buffer.clear();
				appendField(5);
				buffer.tab().append(MEMBERS[m]).tab().append(start).tab().append(stop).tab();
				appendField(8);
				buffer.tab().append(assembly.ordinal).nl();
				for(int i=0; i<writers.length; i++){
					if(assembly.ordinal<panelSizes[i]){writers[i].print(buffer);}
				}
			}
		}
		for(final Assembly assembly:assemblies){
			for(int m=0; m<3; m++){
				require(assembly.observedRaw[m]==assembly.raw[m] && assembly.observedEligible[m]==assembly.eligible[m],
					"Raw truth must retain truncations and conserve saved counts: "+assembly.id+" "+MEMBERS[m]);
			}
		}
	}

	/** Writes ordered source bindings and positive-genome denominators on distinct axes. */
	private void writePanels(){
		try(Output summary=output("panels.tsv")){
			summary.writer.println("panel_genomes\tmember\traw_positive_genomes\tconstruction_positive_genomes\traw_member_loci\tconstruction_originals\tscored_member_positive_genomes");
			for(final int size:panelSizes){
				final long[] raw=new long[3], eligible=new long[3], rawPositive=new long[3], eligiblePositive=new long[3];
				try(Output genomes=output("panel"+size+".genomes.txt"); Output sources=output("panel"+size+".sources.tsv")){
					sources.writer.println("input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff");
					for(int i=0; i<size; i++){
						final Assembly assembly=assemblies.get(i);
						assert(assembly.ordinal==i) : "The same saved prefix order defines genome list and qualified truth";
						genomes.writer.println(assembly.fasta);
						sources.writer.print(buffer.clear().append(i).tab().append(assembly.id).tab().append(assembly.input)
							.tab().append(assembly.fasta).tab().append(assembly.gff).nl());
						for(int m=0; m<3; m++){
							raw[m]+=assembly.raw[m];
							eligible[m]+=assembly.eligible[m];
							rawPositive[m]+=(assembly.raw[m]>0 ? 1 : 0);
							eligiblePositive[m]+=(assembly.eligible[m]>0 ? 1 : 0);
						}
					}
				}
				for(int m=0; m<3; m++){
					summary.writer.print(buffer.clear().append(size).tab().append(MEMBERS[m]).tab().append(rawPositive[m]).tab()
						.append(eligiblePositive[m]).tab().append(raw[m]).tab().append(eligible[m]).append("\tNA_UNION_REFERENCE_ONLY\n"));
				}
			}
		}
	}

	/** Proves the skipped prefix is exactly the already-called population, in order. */
	private void validateCompletedPrefix(){
		int ordinal=0;
		try(Input in=new Input(completedPrefix)){
			RrnaResourceIO.header(in.reader, BINDING_HEADER, completedPrefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(ordinal<chunkStart && p.terms()==5, "Completed prefix must contain exactly chunkstart complete source rows");
				final Assembly assembly=assemblies.get(ordinal);
				require(p.parseInt(0)==ordinal && p.termEquals(assembly.id, 1) && p.parseInt(2)==assembly.input
					&& p.termEquals(assembly.fasta, 3) && p.termEquals(assembly.gff, 4),
					"Completed prefix must match saved order, original input IDs and both source paths at ordinal "+ordinal);
				ordinal++;
			}
		}
		require(ordinal==chunkStart, "Completed prefix ended before the requested skipped population: "+ordinal+" != "+chunkStart);
	}

	/** Writes disjoint suffix lists; ordinals restart per caller, input IDs never change. */
	private void writeChunks(){
		int written=0;
		try(Output summary=output("chunks.tsv")){
			summary.writer.println("chunk\tfirst_ordinal0\tend_ordinal_exclusive\tgenomes\tgenome_list\tsource_bindings");
			for(int start=chunkStart; start<assemblies.size();){
				final int stop=(int)Math.min((long)assemblies.size(), (long)start+chunkSize);
				final String name="chunk"+start;
				final String genomesPath=path(name+".genomes.txt");
				final String sourcesPath=path(name+".sources.tsv");
				try(Output genomes=new Output(genomesPath); Output sources=new Output(sourcesPath)){
					sources.writer.println(BINDING_HEADER);
					for(int i=start; i<stop; i++){
						final Assembly assembly=assemblies.get(i);
						assert(assembly.ordinal==i) : "Disjoint suffix membership must follow the validated full-roster permutation";
						genomes.writer.println(assembly.fasta);
						sources.writer.print(buffer.clear().append(i-start).tab().append(assembly.id).tab().append(assembly.input)
							.tab().append(assembly.fasta).tab().append(assembly.gff).nl());
						written++;
					}
				}
				summary.writer.print(buffer.clear().append(name).tab().append(start).tab().append(stop).tab().append(stop-start)
					.tab().append(genomesPath).tab().append(sourcesPath).nl());
				start=stop;
			}
		}
		require(written==assemblies.size()-chunkStart, "Suffix chunks must conserve every remaining genome, including negatives");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private Assembly assembly(final int inputField, final int idField){
		final int input=p.parseInt(inputField);
		final Assembly assembly=byInput.get(input);
		require(assembly!=null && p.termEquals(assembly.id, idField), "Metadata outside the exact selected input roster: "+input);
		return assembly;
	}
	private int member(final int field){
		for(int i=0; i<MEMBERS.length; i++){if(p.termEquals(MEMBERS[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown producing-CM member: "+p.parseString(field));
	}
	private void appendField(final int field){
		p.appendTerm(buffer, field);
	}
	private Output output(final String name){return new Output(path(name));}
	private String path(final String name){return new File(directory, name).getPath();}
	private static long sum(final int[] counts){
		long total=0;
		for(final int count:counts){total+=count;}
		return total;
	}
	private static int[] parseCounts(final String value, final int expected){
		require(value!=null, "An explicit comma-separated count list is required");
		final LineParser1 parser=new LineParser1(',');
		parser.set(value.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		require(parser.terms()>0 && (expected<0 || parser.terms()==expected), "Unexpected count-list width: "+value);
		final int[] counts=new int[parser.terms()];
		for(int i=0; i<counts.length; i++){counts[i]=parser.parseInt(i);}
		return counts;
	}
	private static void require(final boolean ok, final String why){
		if(!ok){throw new IllegalArgumentException(why);}
	}

	private static final class Assembly{
		Assembly(final String id_, final int input_, final int ordinal_, final int admitted_, final String fasta_, final String source_){
			id=id_;
			input=input_;
			ordinal=ordinal_;
			admitted=admitted_;
			fasta=fasta_;
			source=source_;
		}
		final String id, fasta, source;
		final int input, ordinal, admitted;
		final int[] raw=new int[3], eligible=new int[3], originals=new int[3], observedRaw=new int[3], observedEligible=new int[3];
		String gff;
		int mask;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String orderPath, manifest, countsPath, lociPath, sourcePath, directory, completedPrefix;
	private int expectedGenomes=-1;
	private int chunkStart, chunkSize;
	private int[] expected, rawExpected, panelSizes;
	private final PrintStream outstream;
	private final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
	private final IntObjectMap<Assembly> byInput=new IntObjectMap<Assembly>(Assembly.class);
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String ORDER_HEADER="assembly_id\tinput_id\tchunk\tadmitted_loci\tlinked_fasta\tsource_fasta";
	private static final String BINDING_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	private static final String MANIFEST_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	private static final String COUNTS_HEADER="assembly_id\tinput_id\tcohort\tmember\tnative_events\tnative_loci\tconstruction_loci\ttruncation_only_loci";
	private static final String SOURCE_HEADER="source_id\tmember\tassembly_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tlength\ttaxid\tchunk";
	private static final String LOCUS_HEADER="assembly_id\tinput_id\tcohort\tsource_fasta\tmember\tcontig\tstart\tstop\tstrand\tnative_events\tuntruncated_events\tconstruction_eligible\texact_scored_links";
}
