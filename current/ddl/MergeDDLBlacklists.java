package ddl;

import java.io.PrintStream;

import cardinality.DynamicDemiLog;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import shared.Timer;

/**
 * Merges multiple DDL blacklist FASTA files into one, deduplicating by kmer.
 * Files are loaded in order; each kmer is output only once with the header
 * from the first file that contained it.
 *
 * Usage: MergeDDLBlacklists in=a.fa,b.fa,c.fa out=merged.fa [k=19]
 *
 * @author Brian Bushnell, Noire
 * @date May 24, 2026
 */
public class MergeDDLBlacklists {

	public static void main(String[] args){
		Timer t=new Timer();

		String[] inFiles=null;
		String out=null;
		int k=19;
		boolean overwrite=false;

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in") || a.equals("in1")){
				inFiles=b.split(",");
			}else if(a.equals("out") || a.equals("out1")){
				out=b;
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("ow") || a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}
		}

		assert(inFiles!=null && inFiles.length>0) : "No input files specified.";
		assert(out!=null) : "No output file specified.";

		final LongHashSet seen=new LongHashSet(1024);
		final FileFormat ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		int totalIn=0, totalOut=0, dupes=0;
		for(String path : inFiles){
			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, false, true);
			final ByteFile bf=ByteFile.makeByteFile(ff);
			String header=null;
			StringBuilder sb=new StringBuilder();
			int fileIn=0, fileNew=0;

			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length<1){continue;}
				if(line[0]=='>'){
					if(header!=null && sb.length()>0){
						long packed=DynamicDemiLog.packKmer(sb.toString().getBytes());
						fileIn++;
						if(seen.add(packed)){
							bsw.print(header+"\n");
							bsw.print(sb.toString()+"\n");
							fileNew++;
						}
						sb.setLength(0);
					}
					header=new String(line);
				}else{
					for(byte c : line){sb.append((char)c);}
				}
			}
			if(header!=null && sb.length()>0){
				long packed=DynamicDemiLog.packKmer(sb.toString().getBytes());
				fileIn++;
				if(seen.add(packed)){
					bsw.print(header+"\n");
					bsw.print(sb.toString()+"\n");
					fileNew++;
				}
			}
			bf.close();
			outstream.println(path+": "+fileIn+" in, "+fileNew+" new, "+(fileIn-fileNew)+" duplicates");
			totalIn+=fileIn;
			totalOut+=fileNew;
			dupes+=(fileIn-fileNew);
		}
		bsw.poisonAndWait();

		t.stop();
		outstream.println("Total: "+totalIn+" in, "+totalOut+" unique, "+dupes+" duplicates removed");
		outstream.println("Wrote "+totalOut+" blacklisted kmers to "+out);
		outstream.println("Time: \t"+t);
	}

	private static final PrintStream outstream=System.err;
}
