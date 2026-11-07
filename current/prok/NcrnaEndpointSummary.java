package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntLongHashMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Outputs;
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
				|| radii.length==0 || expected<-1 || expected5<-1 || expected3<-1 || expectedBoth<-1){
			throw new IllegalArgumentException("Read the complete unpaired report; expected counts must be nonnegative or omitted");
		}
		Arrays.sort(radii);
		for(int i=0; i<radii.length; i++){
			if(radii[i]<0 || i>0 && radii[i]==radii[i-1]){throw new IllegalArgumentException("Radii must be distinct nonnegative base counts");}
		}
		reachable=new long[3][radii.length];
		final String[] paths={in, out+".summary.tsv", out+".histogram.tsv", out+".radii.tsv", out+".quantiles.tsv"};
		if(!Tools.testInputFiles(false, true, in) || !Tools.testForDuplicateFiles(true, paths)
				|| !Tools.testOutputFiles(parser.overwrite, false, false, Arrays.copyOfRange(paths, 1, paths.length))){
			throw new IllegalArgumentException("Input/output paths are unreadable, duplicated or not writable");
		}
		input=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		outputs=new FileFormat[4];
		for(int i=0; i<outputs.length; i++){
			outputs[i]=FileFormat.testOutput(paths[i+1], FileFormat.TXT, null, true, parser.overwrite, false, false);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Matched Rows         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates the entire selected population before opening any result writer. */
	public void process(){
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
		checkExpected(expected, count, "matched");
		checkExpected(expected5, Math.max(0, signed[0].get(0)), "exact5");
		checkExpected(expected3, Math.max(0, signed[1].get(0)), "exact3");
		checkExpected(expectedBoth, both, "both-exact");
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
			quantiles.print("family\tend\tmeasure\tmatched\tminimum\tp50\tp90\tp95\tp99\tmaximum\n");
			for(int end=0; end<2; end++){
				final IntLongHashMap absolute=new IntLongHashMap();
				for(final int key:signed[end].toArray()){absolute.increment(Math.abs(key), signed[end].get(key));}
				writeDistribution(histogram, quantiles, signed[end], end, "signed");
				writeDistribution(histogram, quantiles, absolute, end, "absolute");
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
		buffer.clear().append(family).tab().append(label).tab().append(measure).tab().append(count);
		if(count==0){buffer.append("\tNA\tNA\tNA\tNA\tNA\tNA");}
		else{
			buffer.tab().append(keys[0]);
			long cumulative=0;
			int index=0;
			for(final int percent:PERCENTILES){
				final long rank=count/100*percent+(count%100*percent+99)/100;
				while(cumulative<rank){cumulative+=values.get(keys[index++]);}
				buffer.tab().append(keys[index-1]);
			}
			buffer.tab().append(keys[keys.length-1]);
		}
		quantiles.print(buffer.nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Native matched-row input and explicitly declared reporting scope. */
	private final String in;
	private String family, overlap;
	/** Input/output descriptors preserve standard native overwrite/compression policy. */
	private final FileFormat input;
	private final FileFormat[] outputs;
	/** Requested symmetric coordinate radii; counters retain joint reachability. */
	private int[] radii={0, 1, 2, 3, 5, 10, 20, 50, 100};
	private final long[][] reachable;
	/** Optional accepted-grade counts; -1 means unspecified. */
	private long expected=-1, expected5=-1, expected3=-1, expectedBoth=-1;
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
