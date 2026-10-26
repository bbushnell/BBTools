package assemble;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;
import java.util.function.Consumer;

import assemble.HollowTipConnections.Connection;
import assemble.HollowBridgeJoiner.Evidence;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Experimental endpoint-indexed Hollow connection discovery from saved reads.
 * Reports candidate geometry; explicit join=t adds read-verified scaffolding,
 * with experimental sequence filling under fillgaps=t. A unique assembly anchor need not be unique in
 * the genome. Candidate support counts observations, not independent molecules.
 * @author Brian Bushnell, Fischl
 */
public final class HollowBridge {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses the command, runs discovery, and closes redirected diagnostics. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final HollowBridge tool=new HollowBridge(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	/** Uses native BBTools option handling and verifies paths before opening output. */
	public HollowBridge(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 && !split[1].equalsIgnoreCase("null") ? split[1] : null;
			if(a.equals("contigs")){assembly=b;}
			else if(a.equals("extra")){
				if(b!=null){for(String path : b.split(",")){extra.add(path);}}
			}else if(a.equals("k") || a.equals("span")){span=Parse.parseIntKMG(b);}
			else if(a.equals("tipbases")){tipBases=Parse.parseIntKMG(b);}
			else if(a.equals("tipmin")){tipMin=Parse.parseIntKMG(b);}
			else if(a.equals("tipmax")){tipMax=Parse.parseIntKMG(b);}
			else if(a.equals("slide")){slide=Parse.parseIntKMG(b);}
			else if(a.equals("minsupport")){minSupport=Parse.parseIntKMG(b);}
			else if(a.equals("join")){join=Parse.parseBoolean(b);}
			else if(a.equals("fillgaps")){fillGaps=Parse.parseBoolean(b);}
			else if(a.equals("outassembly")){outAssembly=b;}
			else if(a.equals("outjoins")){outJoins=b;}
			else if(a.equals("gapratio")){gapRatio=Parse.parseIntKMG(b);}
			else if(a.equals("flankcheck")){flankCheck=Parse.parseBoolean(b);}
			else if(a.equals("flankbases")){flankBases=Parse.parseIntKMG(b);}
			else if(a.equals("flankmin")){flankMin=Parse.parseIntKMG(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		Parser.processQuality();
		in1=parser.in1;
		in2=parser.in2;
		output=parser.out1;
		overwrite=parser.overwrite;
		threadsIn=parser.threadsIn;
		maxReads=parser.maxReads;
		if(assembly==null || in1==null || output==null || span<=64 || tipBases<32 || minSupport<1){
			throw new IllegalArgumentException("Require contigs=FASTA in=reads out=TSV, k>64, tipbases>=32, minsupport>=1.");
		}
		if(tipMax==-1){tipMax=tipBases-32;}
		if(tipMin<0 || tipMax<tipMin || tipMax>tipBases-32 || slide<0 || (long)tipBases+slide>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Require0<=tipmin<=tipmax<=tipbases-32 and nonnegative slide; tipmax=-1 uses tipbases-32.");
		}
		if(parser.append || parser.out2!=null){throw new IllegalArgumentException("Append and out2 are not supported for a candidate report.");}
		if(flankMin<1 || flankBases<flankMin || flankBases>4096 || (flankCheck && !join)){
			throw new IllegalArgumentException("Require 1<=flankmin<=flankbases<=4096; flankcheck=t requires join=t.");
		}
		if(gapRatio<2 || (join && (outAssembly==null || outJoins==null)) || (!join && (outAssembly!=null || outJoins!=null))){
			throw new IllegalArgumentException("join=t requires outassembly=FASTA and outjoins=TSV; gapratio must be>=2.");
		}
		final ArrayList<String> inputs=new ArrayList<String>();
		inputs.add(in1);
		if(in2!=null){inputs.add(in2);}
		inputs.addAll(extra);
		TadpoleGraph.checkPaths(assembly, inputs, output, outAssembly, outJoins);
		inputs.add(assembly);
		Tools.testInputFiles(false, true, inputs);
		checkInputAliases(inputs);
		if(!Tools.testOutputFiles(overwrite, false, false, output, outAssembly, outJoins)){
			throw new IllegalArgumentException("Output exists or cannot be written: "+output);
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts saved-read evidence, discovers connections, and emits a deterministic report. */
	private void process(final Timer timer){
		assert(span>64 && minSupport>0) : "Discovery requires validated Hollow span and support thresholds.";
		final ArrayList<Contig> contigs=TadpoleGraph.readContigs(assembly);
		final HollowTipConnections finder=new HollowTipConnections(contigs, span, tipBases, minSupport, tipMin, tipMax, slide);
		outstream.println("Assembly contigs: "+contigs.size()+"; unique tip anchors: "+finder.anchors.uniqueCount());
		outstream.println("Counting uses one consumer; t/threads affect native I/O and compression, not table writers.");
		outstream.println("Seed distances: "+tipMin+".."+tipMax+"; sliding radius: "+slide);
		load(in1, in2, finder::addRead);
		for(String path : extra){load(path, null, finder::addRead);}
		final ArrayList<Connection> connections=finder.discover();
		writeReport(contigs, finder, connections);
		int joined=0;
		if(join){
			final HollowBridgeJoiner joiner=new HollowBridgeJoiner(contigs, finder, connections,
					minSupport, gapRatio, flankCheck ? flankBases : 0, flankMin, fillGaps);
			outstream.println(fillGaps ? "Filling accepted joins with read-voted sequence (experimental)."
					: "Scaffolding accepted joins with 10 Ns; complete input contigs retained without overlap trimming.");
			load(in1, in2, joiner::addRead);
			for(String path : extra){load(path, null, joiner::addRead);}
			writeAssembly(joiner.join());
			writeJoins(joiner.evidence());
			joined=joiner.joins();
		}
		int conflicted=0;
		for(Connection connection : connections){if(connection.conflicted()){conflicted++;}}
		timer.stop();
		outstream.println("Observed windows: "+finder.windows()+"; retained: "+finder.retainedWindows());
		outstream.println("Distinct arm pairs: "+finder.distinctPairs()+"; pair array bytes: "+finder.pairArrayBytes());
		outstream.println("Candidates: "+connections.size()+"; conflicted: "+conflicted+"; joins performed: "+joined);
		outstream.println(Tools.timeReadsBasesProcessed(timer, finder.reads(), finder.readBases(), 8));
	}

	/** Streams roots and mates separately; no Hollow window ever spans two reads. */
	private void load(final String first, final String second, final Consumer<byte[]> consumer){
		assert(first!=null && consumer!=null) : "Evidence loading needs a primary file and one mutable count/vote owner.";
		final FileFormat ff1=FileFormat.testInput(first, FileFormat.FASTQ, null, true, true);
		final FileFormat ff2=FileFormat.testInput(second, FileFormat.FASTQ, null, true, true);
		final Streamer reader=StreamerFactory.makeStreamer(ff1, ff2, true, maxReads, false, true, threadsIn);
		reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read : batch.list){
					if(read==null){continue;}
					consumer.accept(read.bases);
					if(read.mate!=null){consumer.accept(read.mate.bases);}
				}
			}
		}finally{reader.close();}
		if(reader.errorState()){throw new RuntimeException("Error reading evidence: "+first);}
	}

	/** Writes complete joined paths and unchanged rejected/isolated contigs through native FASTA output. */
	private void writeAssembly(final ArrayList<Contig> contigs){
		assert(join && outAssembly!=null) : "Assembly output is available only under explicit join=t.";
		final ByteStreamWriter writer=new ByteStreamWriter(outAssembly, overwrite, false, true, FileFormat.FASTA);
		writer.start();
		try{for(Contig contig : contigs){writer.println(contig);}}
		finally{if(writer.poisonAndWait()){throw new RuntimeException("Error writing joined assembly: "+outAssembly);}}
	}

	/** Records every candidate's final disposition, including failed sequence evidence and untouched cycles. */
	private void writeJoins(final ArrayList<Evidence> evidence){
		assert(join && outJoins!=null) : "Join decisions require an explicit output path.";
		final ByteStreamWriter writer=new ByteStreamWriter(outJoins, overwrite, false, true);
		writer.start();
		final ByteBuilder bb=new ByteBuilder();
		try{
			writer.println("#left_end_code\tright_end_code\tgap\tverified_reads\trunner_reads\tflank_rejected\trepeated_pair_windows\tdecision\touter_confirmed\touter_contradictory");
			for(Evidence e : evidence){
				bb.clear();
				bb.append(e.connection.leftEnd).tab().append(e.connection.rightEnd).tab().append(e.connection.gap).tab();
				bb.append(e.support()).tab().append(e.runner()).tab().append(e.flankRejected()).tab().append(e.repeated()).tab().append(e.reason());
				bb.tab().append(e.outerConfirmed()).tab().append(e.outerContradictory());
				writer.println(bb);
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Error writing join decisions: "+outJoins);}}
	}

	/** Writes zero-based input indices/anchor starts and directed physical ends, not inferred genomic coordinates. */
	private void writeReport(final ArrayList<Contig> contigs, final HollowTipConnections finder,
			final ArrayList<Connection> connections){
		assert(contigs!=null && connections!=null) : "Report rows must refer to the same immutable assembly snapshot.";
		final ByteStreamWriter writer=new ByteStreamWriter(output, overwrite, false, true);
		writer.start();
		final ByteBuilder bb=new ByteBuilder();
		try{
			writer.println("#left_index\tleft_name\tleft_end\tright_index\tright_name\tright_end\tgap\tmax_pair_count\tanchor_pairs\tconflicted\tleft_anchor_start\tright_anchor_start\tspan\tslide_checked\tslide_min_count\tslide_failure_offset\tslide_status");
			for(Connection c : connections){
				bb.clear();
				bb.append(c.leftEnd/2).tab();
				appendName(bb, contigs.get(c.leftEnd/2).name);
				bb.tab().append((c.leftEnd&1)==0 ? "left" : "right").tab().append(c.rightEnd/2).tab();
				appendName(bb, contigs.get(c.rightEnd/2).name);
				bb.tab().append((c.rightEnd&1)==0 ? "left" : "right").tab().append(c.gap).tab();
				bb.append(c.support()).tab().append(c.anchorPairs()).tab().append(c.conflicted() ? 1 : 0).tab();
				bb.append(finder.anchors.start(c.leftAnchor())).tab().append(finder.anchors.start(c.rightAnchor())).tab().append(span);
				bb.tab().append(c.slideChecked()).tab().append(c.slideMinCount()).tab().append(c.slideFailureOffset()).tab();
				bb.append(c.slideReason()!=null ? c.slideReason() : slide>0 ? "passed" : "disabled");
				writer.println(bb);
			}
		}finally{
			if(writer.poisonAndWait()){throw new RuntimeException("Error writing candidate report: "+output);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents one physical read file being counted twice through aliases or as assembly evidence. */
	private static void checkInputAliases(final ArrayList<String> paths){
		assert(paths!=null) : "Alias checks require the complete evidence and assembly path list.";
		try{
			final ArrayList<File> files=new ArrayList<File>();
			for(String path : paths){
				if(FileFormat.isStdio(path)){throw new IllegalArgumentException("Use named files for reproducible evidence: "+path);}
				final File file=new File(path).getCanonicalFile();
				for(File previous : files){
					if(previous.equals(file) || Files.isSameFile(previous.toPath(), file.toPath())){
						throw new IllegalArgumentException("Repeated physical input: "+path);
					}
				}
				files.add(file);
			}
		}catch(IOException e){throw new RuntimeException("Could not verify evidence input aliases.", e);}
	}

	/** Retains headers while replacing delimiters that would corrupt the TSV schema. */
	private static void appendName(final ByteBuilder bb, final String name){
		assert(bb!=null) : "Header escaping appends to the report's reusable row buffer.";
		if(name==null){bb.append('.'); return;}
		for(int i=0; i<name.length(); i++){
			final char c=name.charAt(i);
			bb.append(c=='\t' || c=='\r' || c=='\n' ? ' ' : c);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String assembly, in1, in2, output, outAssembly, outJoins;
	private final ArrayList<String> extra=new ArrayList<String>();
	private int span=224, tipBases=192, minSupport=3, threadsIn=-1;
	private int tipMin=0, tipMax=-1, slide=0;
	private int gapRatio=4;
	private int flankBases=128, flankMin=16;
	private boolean flankCheck=false;
	private long maxReads=-1;
	private boolean overwrite=false, join=false, fillGaps=false;
	private PrintStream outstream=System.err;
}
