package prok;

import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Counts RF00174 proposal provenance in saved, final CallGenes GFF rows.
 * Derived from template.A_SampleByteFile. One native byte stream, a reused
 * parser/builder and twelve primitive counters; no per-feature field objects.
 * Every target row survives, including duplicates and unmatched predictions.
 * The original physical line is the join key for BindNcrnaCalls.line_map;
 * contig names alone never establish the owning input. No sequence reads,
 * matching, calling, acceptance replay or biological decisions occur here.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlFinalCallSummary{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, NcrnaDdlFinalCallSummary.class, false);
		try{
			Shared.setThreads(1);
			final Parser parser=new Parser();
			parser.overwrite=false;
			String mode=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("mode")){mode=b;}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown final-call summary parameter: "+arg);}
			}
			require("public".equals(mode) || "hybrid".equals(mode) || "sketch".equals(mode),
				"Declare mode=public,hybrid,sketch to enforce the actual caller arm's metadata contract");
			require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1 && Shared.threads()==1,
				"Supply in=saved_final_GFF out=fresh_prefix without overwrite, append, sampling or parallel readers");
			require(Tools.testInputFiles(false, true, parser.in1), "Saved final CallGenes GFF must be readable");
			Tools.testForDuplicateFiles(true, parser.in1, parser.out1+".calls.tsv", parser.out1+".summary.tsv");
			final NcrnaDdlFinalCallSummary tool=new NcrnaDdlFinalCallSummary(parser.in1, parser.out1, mode);
			tool.process();
			timer.stop();
			pp.outstream.println("Retained RF00174 GFF rows: "+tool.calls+"; source physical lines: "+tool.lines+". Time: "+timer);
			pp.outstream.println("NCRNA_DDL_FINAL_CALL_SUMMARY_COMPLETE_SAVED_GFF_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaDdlFinalCallSummary(final String input_, final String prefix_, final String mode_){
		input=input_;
		prefix=prefix_;
		mode=mode_;
		RrnaResourceIO.fresh(prefix, ".calls.tsv", ".summary.tsv");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		try(Input in=new Input(input); Output out=new Output(prefix+".calls.tsv")){
			out.writer.println("source_physical_line\tmodel\tfirst_proposal_stream\tnative_shared\tnative_support\tseqid\tstart\tstop\tstrand");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				lines++;
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				check(p.terms()==9 && p.termEquals("BBTools", 1), "Expected a nine-field final CallGenes feature");
				if(!p.termEquals("RNA", 2)){continue;}
				attributes(row);
				if(!starts(row, modelStart, modelStop, TARGET)){continue;}
				final int after=modelStart+TARGET.length;
				if(after<modelStop && row[after]!='_' && row[after]!=' ' && row[after]!='\t'){continue;}
				final long start=p.parseLong(3), stop=p.parseLong(4);
				check(start>0 && stop>=start && stop<=Integer.MAX_VALUE && (p.termEquals("+", 6) || p.termEquals("-", 6)),
					"Target calls require positive closed coordinates and explicit strand");
				final int stream=stream(row);
				final int nativeCount=nativeCount(row, stream);
				final int support=nativeCount<0 ? 2 : nativeCount<2 ? 0 : 1;
				counts[stream*3+support]++;
				calls++;
				b.clear().append(lines).tab().append(row, modelStart, modelStop-modelStart)
					.tab().append(STREAMS[stream]).tab();
				if(nativeCount<0){b.append("NA");}
				else{b.append(nativeCount);}
				b.tab().append(SUPPORT[support]).tab();
				field(0).tab().append(start).tab().append(stop).tab();
				field(6).nl();
				out.writer.print(b);
			}
		}
		long total=0;
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.println("first_proposal_stream\tnative_support\tfinal_retained_GFF_rows\tqualification");
			for(int stream=0; stream<STREAMS.length; stream++){
				for(int support=0; support<SUPPORT.length; support++){
					final long count=counts[stream*3+support];
					total=Math.addExact(total, count);
					if(count==0){continue;}
					out.writer.print(b.clear().append(STREAMS[stream]).tab().append(SUPPORT[support]).tab().append(count)
						.append("\tFINAL_CALL_ROWS_NOT_UNIQUE_LOCI_OR_VERIFIER_ACCEPT_EVENTS\n"));
				}
			}
		}
		check(total==calls, "Every retained target row must contribute once; duplicate predictions cannot disappear");
	}

	/** Native Orf attributes use comma/semicolon separators and literal values.
	 * Preserve exact byte spans; a percent sign does not request GFF3 decoding. */
	private void attributes(final byte[] row){
		modelStart=modelStop=streamStart=streamStop=countStart=countStop=-1;
		p.setBounds(8);
		final int limit=p.b();
		for(int from=p.a(); from<limit;){
			int to=from;
			while(to<limit && row[to]!=',' && row[to]!=';'){to++;}
			if(starts(row, from, to, MODEL)){
				check(modelStart<0, "Duplicate model: attribute changes the target identity");
				modelStart=from+MODEL.length;
				modelStop=to;
			}else if(starts(row, from, to, STREAM)){
				check(streamStart<0, "Duplicate ddl_stream attribute is ambiguous");
				streamStart=from+STREAM.length;
				streamStop=to;
			}else if(starts(row, from, to, COUNT)){
				check(countStart<0, "Duplicate ddl_native_shared attribute is ambiguous");
				countStart=from+COUNT.length;
				countStop=to;
			}
			from=to+1;
		}
		check(modelStart>=0 && modelStop>modelStart, "Every generic-RNA row needs its literal model: attribute");
	}

	private int stream(final byte[] row){
		if(mode.equals("public")){
			check(streamStart<0 && countStart<0, "Public calls must remain unclassified; experimental fields contradict mode=public");
			return 3;
		}
		check(streamStart>=0 && countStart>=0, "Experimental target calls require both proposal fields; missing evidence cannot become native support zero");
		for(int stream=0; stream<3; stream++){
			if(equal(row, streamStart, streamStop, STREAM_BYTES[stream])){
				check(mode.equals("sketch") ? stream==2 : stream<2, "Proposal stream contradicts the declared hybrid/sketch caller arm");
				return stream;
			}
		}
		throw failure("Unknown experimental proposal stream");
	}

	private int nativeCount(final byte[] row, final int stream){
		assert(stream>=0 && stream<STREAMS.length) : "Validated proposal streams define the native-support contract";
		if(stream==3){return -1;}
		if(stream==2){
			check(equal(row, countStart, countStop, NA), "Sketch-only native support must be NA, never an invented count");
			return -1;
		}
		check(countStop>countStart && countStop-countStart<=5, "Measured native support must be an unsigned16-bit integer");
		int value=0;
		for(int i=countStart; i<countStop; i++){
			final int digit=row[i]-'0';
			check(digit>=0 && digit<=9, "Hybrid native support requires literal decimal digits");
			value=value*10+digit;
		}
		check(value<=65535, "Measured native support exceeds the caller's unsigned16-bit counter");
		return value;
	}

	private ByteBuilder field(final int field){
		p.setBounds(field);
		return b.append(p.line(), p.a(), p.b()-p.a());
	}
	private static boolean starts(final byte[] row, final int from, final int to, final byte[] prefix){
		assert(from>=0 && from<=to && to<=row.length) : "Attribute comparisons must stay inside the native nine-field row";
		if(to-from<prefix.length){return false;}
		for(int i=0; i<prefix.length; i++){if(row[from+i]!=prefix[i]){return false;}}
		return true;
	}
	private static boolean equal(final byte[] row, final int from, final int to, final byte[] value){
		return to-from==value.length && starts(row, from, to, value);
	}
	private void check(final boolean valid, final String reason){if(!valid){throw failure(reason);}}
	private IllegalArgumentException failure(final String reason){return new IllegalArgumentException(reason+" at "+input+":"+lines);}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static byte[] ascii(final String text){return text.getBytes(StandardCharsets.US_ASCII);}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String input, prefix, mode;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final long[] counts=new long[12];
	private long lines, calls;
	private int modelStart, modelStop, streamStart, streamStop, countStart, countStop;
	private static final byte[] MODEL=ascii("model:"), TARGET=ascii("cobalamin_RF00174"),
		STREAM=ascii("ddl_stream="), COUNT=ascii("ddl_native_shared="), NA=ascii("NA");
	private static final String[] STREAMS={"native", "sketch", "sketch_only", "unclassified"};
	private static final String[] SUPPORT={"BELOW_SHIPPED_FLOOR_2", "AT_LEAST_SHIPPED_FLOOR_2", "NA"};
	private static final byte[][] STREAM_BYTES={ascii("native"), ascii("sketch"), ascii("sketch_only")};
}
