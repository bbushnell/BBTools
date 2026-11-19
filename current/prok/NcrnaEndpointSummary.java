package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntLongHashMap;
import map.ObjectIntMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Outputs;
import prok.NcrnaReferenceGeometry.Target;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Summarizes one family's existing nine-column matched-locus report.
 * Adapted from template.A_SampleByteFile: native byte input, standard a/b CLI
 * parsing and owned native writers. No calling, rematching or truth changes.
 * Positive offsets follow transcript direction. Radius counts describe coordinate
 * reachability, not whether a caller window or neural network can recover an end.
 * Retains only unique locus keys, primitive histograms and radius counters.
 * Optional source bindings exclude exact evaluation assemblies after validating
 * their matched pairs; all retained zero and distant signed errors stay present.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointSummary{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, NcrnaEndpointSummary.class, false);
		try{
			final NcrnaEndpointSummary tool=new NcrnaEndpointSummary(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("Matched loci: "+tool.count);
			pp.outstream.println("Time: "+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses one selected family, declared matching rule and output prefix. */
	public NcrnaEndpointSummary(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("family")){
				family=b;
			}else if(a.equals("bindings")){
				if(b==null || b.isEmpty()){throw new IllegalArgumentException("bindings= requires a source-map path");}
				bindingsPath=b;
			}else if(a.equals("exclude")){
				if(b==null || b.isEmpty()){throw new IllegalArgumentException("exclude= requires a panel-audit path");}
				excludePath=b;
			}else if(a.equals("expectedall")){
				expectedAll=Long.parseLong(b);
			}else if(a.equals("overlap")){
				overlap=b;
			}else if(a.equals("radii")){
				final String[] terms=b.split(",");
				radii=new int[terms.length];
				for(int i=0; i<terms.length; i++){radii[i]=Integer.parseInt(terms[i]);}
			}else if(a.equals("expected")){
				expected=Long.parseLong(b);
			}else if(a.equals("exact5")){
				expected5=Long.parseLong(b);
			}else if(a.equals("exact3")){
				expected3=Long.parseLong(b);
			}else if(a.equals("both")){
				expectedBoth=Long.parseLong(b);
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		in=parser.in1;
		final String out=parser.out1;
		if(in==null || out==null || family==null || family.isEmpty()
				|| !("reciprocal_half".equals(overlap) || "half_shorter".equals(overlap))){
			throw new IllegalArgumentException("Require in=matches.tsv out=prefix family=label overlap=reciprocal_half|half_shorter");
		}
		if(parser.in2!=null || parser.out2!=null || parser.append || parser.maxReads>=0 || parser.samplerate!=1
				|| radii.length==0 || expected<-1 || expected5<-1 || expected3<-1 || expectedBoth<-1 || expectedAll<-1){
			throw new IllegalArgumentException("Read the complete unpaired report; expected counts must be nonnegative or omitted");
		}
		if((bindingsPath==null)!=(excludePath==null) || bindingsPath==null && expectedAll>=0){
			throw new IllegalArgumentException("Source filtering requires bindings= and exclude= together; expectedall= applies before this filter");
		}
		Arrays.sort(radii);
		for(int i=0; i<radii.length; i++){
			if(radii[i]<0 || i>0 && radii[i]==radii[i-1]){throw new IllegalArgumentException("Radii must be distinct nonnegative base counts");}
		}
		reachable=new long[3][radii.length];
		final String[] paths=bindingsPath==null ?
			new String[]{in, out+".summary.tsv", out+".histogram.tsv", out+".radii.tsv", out+".quantiles.tsv"} :
			new String[]{in, out+".summary.tsv", out+".histogram.tsv", out+".radii.tsv", out+".quantiles.tsv", out+".population.tsv"};
		final String[] distinct=Arrays.copyOf(paths, paths.length+(bindingsPath==null ? 0 : 2));
		if(bindingsPath!=null){
			distinct[paths.length]=bindingsPath;
			distinct[paths.length+1]=excludePath;
		}
		if(!Tools.testInputFiles(false, true, in, bindingsPath, excludePath) || !Tools.testForDuplicateFiles(true, distinct)
				|| !Tools.testOutputFiles(parser.overwrite, false, false, Arrays.copyOfRange(paths, 1, paths.length))){
			throw new IllegalArgumentException("Input/output paths are unreadable, duplicated or not writable");
		}
		input=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		outputs=new FileFormat[paths.length-1];
		for(int i=0; i<outputs.length; i++){
			outputs[i]=FileFormat.testOutput(paths[i+1], FileFormat.TXT, null, true, parser.overwrite, false, false);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Matched Rows         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates the entire selected population before opening any result writer. */
	public void process(){
		if(bindingsPath!=null){sources=new SourceFilter(bindingsPath, excludePath);}
		try(Input owner=new Input(input)){
			final ByteFile reader=owner.reader;
			if(!Arrays.equals(reader.nextLine(), HEADER)){throw new IllegalArgumentException("Expected the native nine-column matches header: "+in);}
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				lineNumber++;
				parser.set(line);
				if(parser.terms()!=9){throw bad("Expected nine tab-separated fields");}
				if(parser.termEquals(family, 0)){add();}
			}
		}
		checkExpected(expectedAll, selected, "selected before source exclusion");
		checkExpected(expected, count, "matched");
		checkExpected(expected5, Math.max(0, signed[0].get(0)), "exact5");
		checkExpected(expected3, Math.max(0, signed[1].get(0)), "exact3");
		checkExpected(expectedBoth, both, "both-exact");
		assert(selected==count+excluded) : "Every validated selected-family match is retained or source-excluded exactly once";
		write();
	}

	/** Checks offsets against coordinates and rejects reused truth or call loci. */
	private void add(){
		final boolean minus=parser.termEquals("-", 2);
		if(!minus && !parser.termEquals("+", 2)){throw bad("Strand must be + or -");}
		final int ts=number(3, false), te=number(4, false), cs=number(5, false), ce=number(6, false);
		if(ts<1 || te<ts || cs<1 || ce<cs){throw bad("Coordinates must be nonempty one-based inclusive intervals");}
		final long common=Math.max(0, Math.min(te, ce)-Math.max(ts, cs)+1);
		final long required=overlap.equals("reciprocal_half") ? Math.max(te-ts+1, ce-cs+1) : Math.min(te-ts+1, ce-cs+1);
		if(2*common<required){throw bad("Paired coordinates do not satisfy the declared overlap rule");}
		final int d5=minus ? te-ce : cs-ts, d3=minus ? ts-cs : ce-te;
		if(number(7, true)!=d5 || number(8, true)!=d3){throw bad("Stored offsets disagree with transcript-oriented paired coordinates");}
		if(!truth.add(key(minus, ts, te)) || !calls.add(key(minus, cs, ce))){throw bad("Matched reference or call locus was consumed more than once");}
		if(selected==Long.MAX_VALUE){throw bad("Selected-locus count exceeds long range");}
		selected++;
		if(sources!=null && sources.excludes(parser)){
			excluded++;
			return;
		}
		if(count==Long.MAX_VALUE){throw bad("Matched-locus count exceeds long range");}
		count++;
		signed[0].increment(d5, 1);
		signed[1].increment(d3, 1);
		if(d5==0 && d3==0){both++;}
		final int a5=Math.abs(d5), a3=Math.abs(d3);
		assert(a5>=0 && a3>=0) : "Differences of positive int coordinates cannot equal Integer.MIN_VALUE";
		for(int i=0; i<radii.length; i++){
			if(a5<=radii[i]){reachable[0][i]++;}
			if(a3<=radii[i]){reachable[1][i]++;}
			if(a5<=radii[i] && a3<=radii[i]){reachable[2][i]++;}
		}
	}

	/** Allocates only the two locus keys retained for one-to-one validation. */
	private String key(final boolean minus, final int start, final int stop){
		parser.setBounds(1);
		if(parser.a()==parser.b()){throw bad("Missing contig identifier");}
		return buffer.clear().append(parser.line(), parser.a(), parser.b()-parser.a())
			.tab().append(minus ? '-' : '+').tab().append(start).tab().append(stop).toString();
	}

	/** Strict decimal parsing without per-field Strings or overflow wrapping. */
	private int number(final int term, final boolean signedValue){
		parser.setBounds(term);
		int from=parser.a();
		final boolean negative=signedValue && from<parser.b() && parser.line()[from]=='-';
		if(negative){from++;}
		try{
			final int value=NcrnaNumbers.unsignedInt(parser.line(), from, parser.b());
			return negative ? -value : value;
		}catch(IllegalArgumentException failure){throw bad("Invalid decimal field "+(term+1)+": "+failure.getMessage());}
	}

	private IllegalArgumentException bad(final String message){return new IllegalArgumentException(in+":"+lineNumber+": "+message);}

	/** Optional expectations bind a diagnostic to an already accepted grade. */
	private static void checkExpected(final long expected, final long actual, final String label){
		if(expected>=0 && expected!=actual){throw new IllegalArgumentException(label+" expected "+expected+", observed "+actual);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try(Outputs owners=new Outputs(outputs.length)){
			final ByteStreamWriter summary=owners.open(outputs[0]), histogram=owners.open(outputs[1]);
			final ByteStreamWriter radius=owners.open(outputs[2]), quantiles=owners.open(outputs[3]);
			summary.print("family\toverlap\tmatched\texact5\texact3\tboth_exact\n");
			summary.print(buffer.clear().append(family).tab().append(overlap).tab().append(count).tab()
				.append(Math.max(0, signed[0].get(0))).tab().append(Math.max(0, signed[1].get(0))).tab().append(both).nl());
			radius.print("family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n");
			for(int i=0; i<radii.length; i++){
				radius.print(buffer.clear().append(family).tab().append(overlap).tab().append(radii[i]).tab().append(count)
					.tab().append(reachable[0][i]).tab().append(reachable[1][i]).tab().append(reachable[2][i]).nl());
			}
			histogram.print("family\tend\tmeasure\toffset_bp\tloci\n");
			quantiles.print("#NCRNA_ENDPOINT_QUANTILES_V2\nfamily\tend\tmeasure\tmatched\tminimum\tp50\tp90\tp95\tp99\tmaximum\n");
			for(int end=0; end<2; end++){
				final IntLongHashMap absolute=new IntLongHashMap();
				for(final int key:signed[end].toArray()){absolute.increment(Math.abs(key), signed[end].get(key));}
				writeDistribution(histogram, quantiles, signed[end], end, "signed");
				writeDistribution(histogram, quantiles, absolute, end, "absolute");
			}
			if(sources!=null){
				final ByteStreamWriter population=owners.open(outputs[4]);
				population.print("family\tselected\tretained\texcluded\tbinding_sources\texcluded_assemblies\taudit_panels\taudit_rows\tunmatched_alias_rows\n");
				population.print(buffer.clear().append(family).tab().append(selected).tab().append(count).tab().append(excluded)
					.tab().append(sources.byInput.size()).tab().append(sources.excludedAssemblies)
					.tab().append(sources.panels).tab().append(sources.rows).tab().append(sources.unmatched).nl());
			}
		}
	}

	/** Emits sorted sparse counts and nearest-rank quantiles without expanding counts. */
	private void writeDistribution(final ByteStreamWriter hist, final ByteStreamWriter quantiles,
			final IntLongHashMap values, final int end, final String measure){
		final String label=end==0 ? "5prime" : "3prime";
		final int[] keys=values.toArray();
		Arrays.sort(keys);
		long total=0;
		for(final int key:keys){
			final long n=values.get(key);
			total+=n;
			hist.print(buffer.clear().append(family).tab().append(label).tab().append(measure).tab().append(key).tab().append(n).nl());
		}
		assert(total==count) : "Every matched locus contributes exactly one endpoint offset: histogram="+total+", matched="+count;
		// Quantiles do not exist for an empty population. The summary and radius
		// tables retain its zero counts; do not invent a numeric missing value.
		if(count==0){return;}
		buffer.clear().append(family).tab().append(label).tab().append(measure).tab().append(count);
		buffer.tab().append(keys[0]);
		long cumulative=0;
		int index=0;
		for(final int percent:PERCENTILES){
			final long rank=count/100*percent+(count%100*percent+99)/100;
			while(cumulative<rank){cumulative+=values.get(keys[index++]);}
			buffer.tab().append(keys[index-1]);
		}
		buffer.tab().append(keys[keys.length-1]);
		quantiles.print(buffer.nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------       Source Exclusion       ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact source joins for the optional mode; never infers taxonomy or aliases.
	 * A reusable native byte probe avoids allocating one input-ID string per match.
	 * Source IDs must be printable ASCII, matching the exporter's byte-oriented
	 * length-prefixed namespace. Contig names remain opaque bytes.
	 */
	private static final class SourceFilter{
		SourceFilter(final String bindings, final String exclusions){
			final ArrayList<String> assemblies=new ArrayList<String>();
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(bindings)){
				RrnaResourceIO.header(input.reader, BINDING_HEADER, bindings);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					require(p.terms()==5 && unsigned(p, 0)==assemblies.size(), "Source binding needs dense ordered ordinals and five fields");
					final String assembly=p.parseString(1);
					require(!assembly.isEmpty() && !assembly.equals("."), "A source must bind an exact assembly identity");
					probe.set(p, 2);
					require(probe.to>probe.from && !byInput.contains(probe), "Input source IDs must be nonempty and unambiguous");
					for(int i=probe.from; i<probe.to; i++){
						require(line[i]>=33 && line[i]<=126, "Length-prefixed input IDs must use printable ASCII");
					}
					byInput.put(probe.copy(), assemblies.size());
					assemblies.add(assembly);
				}
			}
			require(!assemblies.isEmpty(), "Source binding must contain the complete nonempty panel");
			final ObjectSet<String> omitted=new ObjectSet<String>(String.class);
			final ObjectIntMap<String> panelRows=new ObjectIntMap<String>(String.class);
			String panel=null;
			try(Input input=new Input(exclusions)){
				RrnaResourceIO.header(input.reader, EXCLUSION_LABEL, exclusions);
				RrnaResourceIO.header(input.reader, EXCLUSION_HEADER, exclusions);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					require(p.terms()==5, "Exclusion audit requires five fields");
					if(panel==null || !p.termEquals(panel, 0)){panel=p.parseString(0);}
					require(panel.startsWith("/") && p.termStartsWith("/", 2), "Exclusion audit retains absolute panel and FASTA paths");
					final int previous=panelRows.get(panel), ordinal=unsigned(p, 1);
					require(ordinal==(previous<0 ? 1 : previous+1), "Exclusion panel row indices must be dense and ordered");
					panelRows.put(panel, ordinal);
					rows++;
					if(p.termEquals("EXCLUDED_ASSEMBLY", 4)){
						final String assembly=p.parseString(3);
						require(!assembly.isEmpty() && !assembly.equals("."), "Excluded rows require exact assembly identities");
						omitted.add(assembly);
					}else{
						require(p.termEquals("UNMATCHED_MANIFEST_PATH_ALIAS", 4) && p.termEquals(".", 3),
							"Unknown exclusion status or contradictory unmatched assembly");
						unmatched++;// Report explicitly; the family worker must bind the accepted unmatched count.
					}
				}
			}
			panels=panelRows.size();
			excludedAssemblies=omitted.size();
			require(panels>0, "Exclusion audit must contain at least one evaluation panel");
			excluded=new boolean[assemblies.size()];
			for(int i=0; i<excluded.length; i++){excluded[i]=omitted.contains(assemblies.get(i));}
		}

		/** Reads exactly length:input_id/contig; length is never an input ordinal. */
		boolean excludes(final LineParser1 p){
			p.setBounds(1);
			final byte[] line=p.line();
			final int from=p.a(), stop=p.b();
			int colon=from;
			while(colon<stop && line[colon]>='0' && line[colon]<='9'){colon++;}
			require(colon>from && colon<stop && line[colon]==':', "Expected length-prefixed matched source ID");
			final int length=NcrnaNumbers.unsignedInt(line, from, colon);
			final long end=(long)colon+1+length;
			require(length>0 && end<stop-1 && line[(int)end]=='/', "Declared source-ID length must end at the contig separator");
			probe.set(line, colon+1, (int)end);
			final int source=byInput.get(probe);
			require(source>=0, "Matched input ID is outside the declared source binding");
			return excluded[source];
		}

		private static int unsigned(final LineParser1 p, final int field){
			p.setBounds(field);
			return NcrnaNumbers.unsignedInt(p.line(), p.a(), p.b());
		}
		private static void require(final boolean valid, final String reason){
			if(!valid){throw new IllegalArgumentException(reason);}
		}
		private final ObjectIntMap<Target> byInput=new ObjectIntMap<Target>(Target.class);
		private final Target probe=new Target();
		private final boolean[] excluded;
		private int excludedAssemblies, panels;
		private long rows, unmatched;
		private static final String BINDING_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
		private static final String EXCLUSION_LABEL="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result";
		private static final String EXCLUSION_HEADER="panel\trow1\tfasta\tassembly\tstatus";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Native matched-row input and explicitly declared reporting scope. */
	private final String in;
	private String family, overlap;
	private String bindingsPath, excludePath;
	private SourceFilter sources;
	/** Input/output descriptors preserve standard native overwrite/compression policy. */
	private final FileFormat input;
	private final FileFormat[] outputs;
	/** Requested symmetric coordinate radii; counters retain joint reachability. */
	private int[] radii={0, 1, 2, 3, 5, 10, 20, 50, 100};
	private final long[][] reachable;
	/** Optional accepted-grade counts; -1 means unspecified. */
	private long expected=-1, expected5=-1, expected3=-1, expectedBoth=-1;
	private long expectedAll=-1, selected, excluded;
	private long count, both, lineNumber=1;
	/** Reused parsing/output scratch; retained sets enforce one-to-one locus use. */
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final IntLongHashMap[] signed={new IntLongHashMap(), new IntLongHashMap()};
	private final ObjectSet<String> truth=new ObjectSet<String>(String.class);
	private final ObjectSet<String> calls=new ObjectSet<String>(String.class);
	private static final int[] PERCENTILES={50, 90, 95, 99};
	private static final byte[] HEADER=("#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3")
		.getBytes(StandardCharsets.UTF_8);
}
