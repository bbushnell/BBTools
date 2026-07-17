package synth;

import java.io.PrintStream;

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
 * Fuses contigs sharing a taxonomy ID (from "tid|TAXID|..." headers) into N-padded genome-scale
 * sequences, emitting only fused sequences at least 'floor' bp long.  This turns a combined RefSeq
 * clade file (mostly short rRNA markers plus some assemblies) into assembly-attempt sequences for
 * CONTIG-classification benchmarking, dropping marker-only taxa.
 *
 * REQUIRES taxonomically-sorted input (contigs of the same taxid adjacent), e.g. from
 * sortbyname.sh taxa=t, so it can stream without holding the whole file in memory.
 *
 * Usage: fusebytaxa.sh in=sorted.fna.gz out=genomes.fna.gz floor=1m pad=10 maxlen=50m
 *   floor=1m   Only emit fused sequences at least this long (0 to emit all, e.g. viruses).
 *   pad=10     Ns inserted between fused contigs.
 *   maxlen=50m Split a taxon's fused sequence into chunks no longer than this.
 *
 * @author Noire
 * @date July 16, 2026
 */
public class FuseByTaxa {

	public static void main(String[] args){
		Timer t=new Timer();
		FuseByTaxa x=new FuseByTaxa(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public FuseByTaxa(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("floor") || a.equals("minlen")){
				floor=Parse.parseKMG(b);
			}else if(a.equals("pad") || a.equals("npad")){
				pad=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxlen=Parse.parseKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//standard flags
			}else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		Parser.processQuality();
		overwrite=parser.overwrite;
		in1=parser.in1;
		out1=parser.out1;
		assert(in1!=null) : "No input file specified.";
		assert(maxlen>=floor) : "maxlen ("+maxlen+") must be >= floor ("+floor+").";

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	void process(Timer t){
		//ordered=true and threadsIn=1 keep contigs in input order, which the fuse-adjacent logic needs.
		Streamer st=StreamerFactory.makeStreamer(ffin1, null, true, -1L, false, true, 1);
		st.start();
		final ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1));
		if(bsw!=null){bsw.start();}

		final ByteBuilder bb=new ByteBuilder();   //current fused chunk (bases only)
		String curTid=null;
		int chunk=0;
		boolean emittedThisTaxon=false;

		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				readsIn++;
				basesIn+=r.length();
				final String tid=parseTid(r.id);
				if(curTid==null){curTid=tid;}
				if(!tid.equals(curTid)){
					//New taxon: flush the current chunk, then account for whether the taxon was kept.
					if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
					if(emittedThisTaxon){taxaKept++;}else{taxaDropped++;}
					bb.clear(); curTid=tid; chunk=0; emittedThisTaxon=false;
				}else if(bb.length()>0 && bb.length()+pad+r.length()>maxlen){
					//Same taxon but the chunk is full: flush and continue a new chunk.
					if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
					bb.clear(); chunk++;
				}
				if(bb.length()>0){for(int k=0; k<pad; k++){bb.append((byte)'N');}}
				bb.append(r.bases);
			}
		}
		//Flush the final taxon.
		if(curTid!=null){
			if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
			if(emittedThisTaxon){taxaKept++;}else{taxaDropped++;}
		}

		st.close();
		if(bsw!=null){bsw.poisonAndWait();}
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println("Taxa kept:    \t"+taxaKept);
		outstream.println("Taxa dropped: \t"+taxaDropped+" (below floor="+floor+")");
		outstream.println("Sequences out:\t"+seqsOut+"  ("+basesOut+" bp)");
	}

	/** Writes the fused chunk as one FASTA record if it meets the floor; returns true if written. */
	private boolean emit(ByteStreamWriter bsw, String tid, int chunk, ByteBuilder bb){
		if(bb.length()<floor || bb.length()==0){return false;}
		seqsOut++;
		basesOut+=bb.length();
		if(bsw!=null){
			header.clear();
			header.append('>').append("tid|").append(tid).append("|fused").append(chunk)
				.append(" len=").append(bb.length()).nl();
			bsw.print(header);
			bsw.print(bb);
			bsw.print(NEWLINE);
		}
		return true;
	}

	/** Extracts the taxid from a "tid|TAXID|..." header; returns "NA" if absent. */
	private static String parseTid(String id){
		if(id==null){return "NA";}
		int a=id.indexOf('|');
		if(a<0){return "NA";}
		int b=id.indexOf('|', a+1);
		return b<0 ? id.substring(a+1) : id.substring(a+1, b);
	}

	/*--------------------------------------------------------------*/

	private String in1=null;
	private String out1=null;
	private FileFormat ffin1;
	private FileFormat ffout1;

	private long floor=1000000;   //minimum fused length to emit
	private int pad=10;           //Ns between fused contigs
	private long maxlen=50000000; //split fused sequences longer than this

	private long readsIn=0, basesIn=0;
	private long seqsOut=0, basesOut=0;
	private long taxaKept=0, taxaDropped=0;

	private final ByteBuilder header=new ByteBuilder(64);
	private static final byte[] NEWLINE={(byte)'\n'};

	private boolean overwrite=false;
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
}
