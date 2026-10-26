package test.assemble;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Evaluates an experimental start veto on saved right-extensions, without rerunning them.
 * A rejected record is truncated to its original seed length; no base is corrected.
 * Truth is deliberately absent from the gate and must be audited independently afterward.
 * @author Brian Bushnell, Fischl */
public final class HollowSpanGate {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses standard config expansion and closes redirected diagnostics. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowSpanGate.class, false);
		try{new HollowSpanGate(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires three ordered terminal reports and fresh sequence/report outputs. */
	private HollowSpanGate(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("minus")){minus=b;}
			else if(a.equals("nominal")){nominal=b;}
			else if(a.equals("plus")){plus=b;}
			else if(a.equals("details")){details=b;}
			else if(a.equals("k")){k=Parse.parseIntKMG(b);}
			else if(a.equals("ratio")){ratio=Parse.parseIntKMG(b);}
			else if(a.equals("minsupport")){minSupport=Parse.parseIntKMG(b);}
			else if(a.equals("test")){test=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		require(k>=67 && ratio>=1 && minSupport>=1, "Require k>=67, ratio>=1 and minsupport>=1.");
		require(!parser.overwrite && !parser.append, "Guard outputs must be new, never overwritten or appended.");
		if(test){return;}
		require(input!=null && output!=null && details!=null && minus!=null && nominal!=null && plus!=null,
				"Specify in=extensions, out=guarded.fa, details=flags.tsv, minus=, nominal= and plus= terminal reports.");
		require(Tools.testInputFiles(false, true, input, minus, nominal, plus), "Missing guard input.");
		require(Tools.testOutputFiles(false, false, false, output, details) && !output.equals(details), "Guard outputs must be distinct and new.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Joins by both strict stream order and ID, checking endpoint identity before any decision. */
	private void process(){
		if(test){selfTest(); return;}
		final Streamer stream=StreamerFactory.makeStreamer(
				FileFormat.testInput(input, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
		final Writer writer=WriterFactory.makeWriter(
				FileFormat.testOutput(output, FileFormat.FASTA, null, true, false, false, true), null, 1, null, false);
		final HashSet<String> seen=new HashSet<String>();
		final ByteBuilder row=new ByteBuilder();
		long reads=0, flagged=0, withheld=0, removed=0;
		boolean started=false, error=false;
		try(Profile low=new Profile(minus); Profile mid=new Profile(nominal); Profile high=new Profile(plus);
				Report report=new Report(details)){
			report.writer.print("id\tseed_length\tminus_depth\tnominal_depth\tplus_depth\toriginal_added\n");
			stream.start(); writer.start(); started=true;
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read read : batch){
					require(read.mate==null && seen.add(read.id), "Unknown pairing or repeated extension ID: "+read.id);
					final int a=low.consume(read, k-1, k), b=mid.consume(read, k, k), c=high.consume(read, k+1, k);
					require(low.seedLength==mid.seedLength && high.seedLength==mid.seedLength, "Profile seed lengths differ: "+read.id);
					reads++;
					if(reject(a, b, c, ratio, minSupport)){
						final int added=read.length()-mid.seedLength;
						flagged++;
						if(added>0){withheld++; removed+=added;}
						row.clear().append(read.id).tab().append(mid.seedLength).tab().append(a).tab().append(b)
								.tab().append(c).tab().append(added).nl();
						report.writer.print(row);
						if(added>0){
							read.bases=Arrays.copyOf(read.bases, mid.seedLength);
							if(read.quality!=null){read.quality=Arrays.copyOf(read.quality, mid.seedLength);}
						}
					}
				}
				writer.addReads(batch);
			}
			low.requireEnd(); mid.requireEnd(); high.requireEnd();
		}finally{
			if(started){error|=writer.poisonAndWait(); error|=ReadWrite.closeStream(stream);}
		}
		require(!error, "Guard sequence I/O failed.");
		System.err.println("HOLLOW_SPAN_GATE_PASS reads="+reads+" flagged="+flagged+" withheld_extensions="+withheld+" removed_bases="+removed);
	}

	/** An exploratory abstention rule, not a next-base choice or a claim of an indel. */
	private static boolean reject(final int minus, final int nominal, final int plus, final int ratio, final int floor){
		assert(minus>=0 && nominal>=0 && plus>=0 && ratio>=1 && floor>=1) : "Depth reports and CLI validation supply nonnegative counts and positive thresholds.";
		final int alternate=Math.max(minus, plus);
		return alternate>=floor && alternate>=(long)ratio*nominal;
	}

	/** Exercises both directions, equality, minimum support and multiplication overflow. */
	private static void selfTest(){
		require(reject(451, 9, 0, 4, 3), "Known insertion-spacing pilot must be flagged.");
		require(!reject(0, 7, 1, 4, 3), "Substitution positive must remain eligible.");
		require(!reject(11, 3, 0, 4, 3), "Below-ratio alternative must not pass.");
		require(reject(12, 3, 0, 4, 3), "Ratio equality is inclusive.");
		require(reject(0, 3, 12, 4, 3), "Longer-spacing alternative must also be tested.");
		require(!reject(2, 0, 0, 4, 3), "Minimum support protects zero-depth cases.");
		require(!reject(Integer.MAX_VALUE, Integer.MAX_VALUE, 0, 4, 3), "Integer multiplication must not overflow.");
		require(!reject(0, 0, 0, 4, 3), "Absent evidence is not a veto.");
		System.out.println("HOLLOW_SPAN_GATE_TEST_PASS checks=8");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** One native tabular stream; only its current row remains live. */
	private static final class Profile implements AutoCloseable {
		/** Validates schema before the record stream starts, closing on constructor failure. */
		Profile(final String path){
			file=ByteFile.makeByteFile(path, true);
			boolean valid=false;
			try{
				final byte[] header=file.nextLine();
				valid=header!=null && new String(header, java.nio.charset.StandardCharsets.US_ASCII).equals(SCHEMA);
				require(valid, "Unexpected terminal profile header: "+path);
			}finally{if(!valid){file.close();}}
		}

		/** Exact order, ID, spans and both anchors must match the original extension prefix. */
		int consume(final Read read, final int tableK, final int span){
			final byte[] line=file.nextLine();
			require(line!=null, "Terminal profile ended before extensions: "+read.id);
			lp.set(line);
			require(lp.terms()==7 && lp.termEquals(read.id, 0), "Profile ID/order mismatch: "+read.id);
			seedLength=lp.parseInt(1);
			require(seedLength>=span && seedLength<=read.length() && lp.parseInt(2)==tableK && lp.parseInt(3)==span,
					"Seed length/table span/query span mismatch: "+read.id);
			anchor(read.bases, seedLength-span, 4);
			anchor(read.bases, seedLength-32, 5);
			final int depth=lp.parseInt(6);
			require(depth>=0, "Negative terminal depth.");
			return depth;
		}

		/** Compares native parser byte ranges without allocating anchor Strings or arrays. */
		private void anchor(final byte[] bases, final int offset, final int field){
			assert(offset>=0 && offset<=bases.length-32) : "consume validates seed bounds before comparing anchors.";
			require(lp.setBounds(field)==32, "Terminal anchors must contain32 bases.");
			for(int i=0; i<32; i++){require(lp.parseByteFromCurrentField(i)==bases[offset+i], "Terminal anchor disagrees with saved seed prefix.");}
		}

		/** Extra rows are fatal rather than silently left unused. */
		void requireEnd(){require(file.nextLine()==null, "Terminal profile has extra rows.");}
		@Override
		public void close(){require(!file.close(), "Terminal profile read failed.");}
		private final ByteFile file;
		private final LineParser1 lp=new LineParser1('\t');
		int seedLength;
	}

	/** Guarantees terminal delivery to the report writer even when a record check fails. */
	private static final class Report implements AutoCloseable {
		/** Opens one fresh asynchronous report owned by this resource. */
		Report(final String path){
			writer=new ByteStreamWriter(path, false, false, true);
			writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Guard report write failed.");}
		final ByteStreamWriter writer;
	}

	/** Keeps scientific input checks active independently of -ea. */
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, output, details, minus, nominal, plus;
	private int k=128, ratio=4, minSupport=3;
	private boolean test=false;
	private static final String SCHEMA="id\tlength\tk\tquery_span\tleft_anchor\tright_anchor\tdepth";
}
