package dna;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextStreamWriter;
import parse.PreParser;
import shared.Timer;
import shared.Tools;

/**
 * Converts chromosome array files to FASTA format.
 * Supports both whole chromosome output and contig splitting based on N-run length.
 * Can process individual chromosome files or entire directories of chromosome files.
 * Directory input retains filesystem enumeration order. Whole-chromosome output
 * uses chromosome numbers as headers; split output numbers contigs consecutively.
 * Helpers taking a writer borrow it and do not close it. The filename entry points
 * own their writer and check output errors after joining on normal completion.
 * All entry points reject nonpositive line widths before opening or writing FASTA,
 * including empty input. CLI minus signs on ASCII integer tokens survive config
 * expansion and option filtering.
 *
 * @author Brian Bushnell
 * @date Jul 26, 2012
 */
public class ChromToFasta{

	/** Converts serialized chromosome arrays using the legacy positional interface.
	 * @param args Input file list/directory/build=N, output path, positive line width,
	 * and optional N-run split threshold; a nonpositive threshold emits whole chromosomes
	 * @throws IllegalArgumentException If the line width is nonpositive, before opening FASTA output
	 * @throws RuntimeException If input loading or checked output completion fails */
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			//DNA-050: retain signed positional integers after config expansion and global-option removal.
			PreParser pp=new PreParser(args, System.err, new Object() { }.getClass().getEnclosingClass(), false, true, true, true);
			args=pp.args;
			//TODO: Probable bug (DNA-052) - outstream=file opens an owned status stream
			//which this main discards without closing; process exit masks the leak in standalone use.
			//outstream=pp.outstream;
		}

		Timer t=new Timer();

		if(args[0].contains("=") && (args[0].startsWith("build") || args[0].startsWith("genome"))){
			int build=Integer.parseInt(args[0].split("=")[1]);
//			Data.setGenome(build);
//			String s="", comma="";
//			for(int i=1; i<Data.numChroms; i++){
//				s=s+comma+Data.chromFname(i, build);
//				comma=",";
//			}
//			args[0]=s;
			args[0]=Data.chromFname(1, build);
			args[0]=args[0].substring(0, args[0].lastIndexOf('/'));
		}

		String[] chromfiles=args[0].split(",");

		if(chromfiles.length==1){
			File f=new File(chromfiles[0]);
			if(f.isDirectory()){
				ArrayList<String> list=new ArrayList<String>(4);
				for(File f2 : f.listFiles()){
					if(!f2.isDirectory() && f2.isFile()){
						String s=f2.getAbsolutePath();
						if(s.endsWith(".chrom") || s.endsWith(".chromC") || s.contains(".chrom.") || s.contains(".chromC.")){
							list.add(s);
						}
					}
				}
				chromfiles=list.toArray(new String[list.size()]);
			}
		}

		String outfile=args[1];
		int blocklen=Integer.parseInt(args[2]);
		if(blocklen<=0){throw new IllegalArgumentException("FASTA line width must be positive: "+blocklen);}
		int trigger=(args.length>3 ? Integer.parseInt(args[3]) : 0);

		//TODO: Probable bug (DNA-020) - input/formatting exceptions after start skip
		//writer shutdown. Preserve this path until a bounded exceptional-cleanup test exists.
		TextStreamWriter tsw=new TextStreamWriter(outfile, true, false, false);
		tsw.start();

		if(trigger<=0){ //Write normally
			for(int i=0; i<chromfiles.length; i++){
				ChromosomeArray cha=ChromosomeArray.read(chromfiles[i]);
				writeChrom(cha, tsw, blocklen);
			}
		}else{ //Break into contigs
			int contig=1;
			for(int i=0; i<chromfiles.length; i++){
				ChromosomeArray cha=ChromosomeArray.read(chromfiles[i]);
				contig=writeContigs(cha, contig, trigger, blocklen, tsw);
			}
		}

		if(tsw.poisonAndWait()){
			throw new RuntimeException("Failed writing FASTA output: "+outfile);
		}

		t.stop();
		System.err.println("Time:\t"+t);
	}

	/** Splits the inclusive stored interval at uppercase N runs of at least trigger bases.
	 * Omits leading/trailing Ns; shorter internal N runs remain in the sequence.
	 * @param cha Source chromosome with inclusive minIndex/maxIndex bounds
	 * @param contig Identifier for the first emitted contig
	 * @param trigger Positive N-run threshold
	 * @param fastaBlocklen Positive output line width
	 * @param tsw Started caller-owned writer; remains open on return
	 * @return Identifier available for the next contig, unchanged if none was emitted
	 * @throws IllegalArgumentException If the line width is nonpositive, before any output */
	public static int writeContigs(ChromosomeArray cha, int contig, int trigger, int fastaBlocklen, TextStreamWriter tsw){
		if(fastaBlocklen<=0){throw new IllegalArgumentException("FASTA line width must be positive: "+fastaBlocklen);}

		StringBuilder sb=new StringBuilder(4000);

		//Note: ns counts consecutive Ns and is NOT reset after an emit — safe because Ns are only
		//appended to sb while sb.length()>0, so trailing Ns of a just-emitted (now-empty) contig are
		//skipped and never mis-trimmed. sb always begins with a non-N base, so the setLength(len-ns)
		//trims exactly the trailing Ns and can never underflow.
		int ns=0;

		for(int aloc=cha.minIndex; aloc<=cha.maxIndex; aloc++){
			byte b=cha.get(aloc);
			if(b=='N'){
				ns++;
				if(sb.length()>0){
					sb.append('N');
					if(ns==trigger){
						sb.setLength(sb.length()-ns);
						tsw.print(">"+contig+"\n");
						contig++;
						writeContig(sb, tsw, fastaBlocklen);
						sb.setLength(0);
					}
				}
			}else{
				sb.append((char)b);
				ns=0;
			}
		}

		if(sb.length()>0){
			sb.setLength(sb.length()-ns);
			tsw.print(">"+contig+"\n");
			contig++;
			writeContig(sb, tsw, fastaBlocklen);
			sb.setLength(0);
		}

		return contig;
	}

	/** Queues wrapped sequence lines without a header; does not mutate the supplied builder.
	 * @param sb Contig sequence; an empty builder emits no lines
	 * @param tsw Started caller-owned writer, left open
	 * @param blocklen Positive line width
	 * @throws IllegalArgumentException If the line width is nonpositive, before any output */
	public static void writeContig(StringBuilder sb, TextStreamWriter tsw, int blocklen){
		if(blocklen<=0){throw new IllegalArgumentException("FASTA line width must be positive: "+blocklen);}
		//TODO: Probable bug (DNA-051) - positive widths can still overflow i+blocklen
		//or the loop increment at large coordinates; separate bounded arithmetic repair needed.
		for(int i=0; i<sb.length(); i+=blocklen){
			int max=Tools.min(i+blocklen, sb.length());
			tsw.println(sb.substring(i, max));
		}
	}

	/** Writes one chromosome and waits for output closure; the delegate emits its single header.
	 * Overwrites the destination. Shutdown is not guaranteed if formatting itself throws.
	 * @param cha Source chromosome
	 * @param fname Output FASTA path
	 * @param blocklen Positive output line width
	 * @throws IllegalArgumentException If the line width is nonpositive, before opening the destination
	 * @throws RuntimeException If writing or closing the FASTA output fails */
	public static void writeChrom(ChromosomeArray cha, String fname, int blocklen){
		if(blocklen<=0){throw new IllegalArgumentException("FASTA line width must be positive: "+blocklen);}
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		//TODO: Probable bug (DNA-020) - a delegate exception bypasses owned-writer shutdown.
		writeChrom(cha, tsw, blocklen);
		if(tsw.poisonAndWait()){
			throw new RuntimeException("Failed writing FASTA output: "+fname);
		}
	}

	/** Queues one chromosome header and bases from coordinate zero through maxIndex.
	 * Coordinates before minIndex are emitted as N via getString; does not close the writer.
	 * @param cha Source chromosome with inclusive maximum coordinate
	 * @param tsw Started caller-owned writer
	 * @param blocklen Positive output line width
	 * @throws IllegalArgumentException If the line width is nonpositive, before any output */
	public static void writeChrom(ChromosomeArray cha, TextStreamWriter tsw, int blocklen){
		if(blocklen<=0){throw new IllegalArgumentException("FASTA line width must be positive: "+blocklen);}
		tsw.println(">"+cha.chromosome);
		//TODO: Probable bug (DNA-051) - i+blocklen-1 and i+=blocklen can overflow
		//with positive widths at large coordinates; width validation does not repair that arithmetic.
		for(int i=0; i<=cha.maxIndex; i+=blocklen){
			int max=Tools.min(i+blocklen-1, cha.maxIndex);
			tsw.println(cha.getString(i, max));
		}
	}

}
