package prot;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntLongHashMap;
import tax.TaxTree;

/**
 * Builds sizemap.tsv (tid&lt;TAB&gt;bp) for the MAG-QC foundation corpus: total
 * genome length (all contigs summed) per tid, the completeness denominator
 * required by MagQCVectorMaker.
 *
 * Scans one or more genome directories for *.fna.gz, parses tid from each
 * filename (tid_&lt;N&gt;_..., via {@link TaxTree#parseTaxID(String)} -- the
 * same convention CacheBuilder uses for its domain-set filenames), and sums
 * sequence bp per tid across ALL contributing files. A tid backed by more than
 * one file is not an error: two files under one tid are two parts of the same
 * organism's genome (confirmed correct by UMP45, 2026-08-17) and their bp is
 * summed, never overwritten. Every such case is reported by name (tid + the
 * contributing filenames), never silently merged.
 *
 * Usage: sizemapbuilder.sh indir=bacteria5 indir2=archaea4 out=sizemap.tsv
 *
 * @author Eru
 */
public class SizeMapBuilder {

	public static void main(String[] args){
		Timer t=new Timer();
		SizeMapBuilder x=new SizeMapBuilder(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public SizeMapBuilder(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("indir") || a.equals("in") || a.equals("indir2") || a.equals("in2")
					|| a.equals("indir3") || a.equals("in3")){
				dirs.add(b);
			}else if(a.equals("out")){outFile=b;}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(!dirs.isEmpty()) : "At least one indir= is required.";
		assert(outFile!=null) : "out= is required.";
	}

	void process(Timer t){
		final IntLongHashMap bpMap=new IntLongHashMap(1<<16);
		final IntHashMap fileCount=new IntHashMap(1<<16);
		final HashMap<Integer, ArrayList<String>> filesByTid=new HashMap<Integer, ArrayList<String>>();
		long filesScanned=0;

		for(String dir : dirs){
			File d=new File(dir);
			File[] files=d.listFiles();
			assert(files!=null) : "Could not list directory: "+dir;
			ArrayList<File> fnaFiles=new ArrayList<File>();
			for(File f : files){
				if(f.getName().endsWith(".fna.gz")){fnaFiles.add(f);}
			}
			assert(!fnaFiles.isEmpty()) : "No .fna.gz files found in "+dir;
			outstream.println("Scanning "+fnaFiles.size()+" .fna.gz files in "+dir+"...");

			for(File f : fnaFiles){
				final String name=f.getName();
				final int tid=TaxTree.parseTaxID(name);
				assert(tid>0) : "Could not parse tid from filename: "+name;

				final long bp=scanBp(f.getPath());
				assert(bp>0) : "Zero-length sequence content in "+f.getPath();

				if(!bpMap.put(tid, bp)){
					// Key already present: IntLongHashMap.put() does not overwrite (returns false),
					// so accumulate explicitly -- remove then put the summed value.
					final long old=bpMap.get(tid);
					bpMap.remove(tid);
					bpMap.put(tid, old+bp);
				}
				final int fc=fileCount.get(tid);
				fileCount.put(tid, (fc<0 ? 1 : fc+1));

				ArrayList<String> names=filesByTid.get(tid);
				if(names==null){names=new ArrayList<String>(2); filesByTid.put(tid, names);}
				names.add(name);

				filesScanned++;
				if(filesScanned%2000==0){
					outstream.println("  "+filesScanned+" files scanned, "+bpMap.size()+" distinct tids");
				}
			}
		}

		// Every multi-file tid must be named, not silently merged (UMP45's requirement).
		int multiFileTids=0;
		for(HashMap.Entry<Integer, ArrayList<String>> e : filesByTid.entrySet()){
			if(e.getValue().size()>1){
				multiFileTids++;
				outstream.println("  tid "+e.getKey()+" backed by "+e.getValue().size()
					+" files (bp summed): "+e.getValue());
			}
		}
		outstream.println("Distinct tids with >1 contributing file: "+multiFileTids);

		outstream.println("Writing "+bpMap.size()+" tid->bp rows to "+outFile+"...");
		final FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		try{
			final ByteBuilder bb=new ByteBuilder(1<<12);
			final int[] keys=bpMap.keys();
			final int invalid=bpMap.invalid();
			int n=0;
			for(int i=0; i<keys.length; i++){if(keys[i]!=invalid){n++;}}
			final int[] sortedKeys=new int[n];
			int ki=0;
			for(int i=0; i<keys.length; i++){if(keys[i]!=invalid){sortedKeys[ki++]=keys[i];}}
			java.util.Arrays.sort(sortedKeys);
			for(int tid : sortedKeys){
				final long bp=bpMap.get(tid);
				assert(bp>0) : "Non-positive summed bp for tid "+tid;
				bb.append(tid).append('\t').append(bp).nl();
				bsw.print(bb); bb.clear();
			}
		}finally{
			bsw.poisonAndWait();
		}

		t.stop();
		outstream.println("Files scanned: "+filesScanned+"   distinct tids: "+bpMap.size());
		outstream.println("Time: \t"+t);
	}

	/** Sums sequence bp (all non-header bytes, including any embedded N) in a gzipped FASTA file. */
	static long scanBp(String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		long bp=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]=='>'){continue;}
			bp+=line.length;
		}
		bf.close();
		return bp;
	}

	private ArrayList<String> dirs=new ArrayList<String>();
	private String outFile=null;
	private boolean overwrite=true;
	private PrintStream outstream=System.err;
}
