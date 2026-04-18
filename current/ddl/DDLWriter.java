package ddl;

import java.io.File;

import cardinality.DynamicDemiLog;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import parse.Parse;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Builds DynamicDemiLog sketches from FASTA/FASTQ files and writes them
 * as A48-encoded TSV. Supports per-file mode; per-sequence and per-id
 * modes are planned.
 *
 * Usage: DDLWriter in=a.fa,b.fa out=ddls.tsv.gz [k=31] [buckets=2048]
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLWriter {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		DDLWriter dw=new DDLWriter(args);
		dw.process(t);
	}

	public DDLWriter(String[] args){
		Shared.capBufferLen(200);
		Shared.capBuffers(8);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("buckets")){
				buckets=Integer.parseInt(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("mode")){
				if(b.equalsIgnoreCase("perfile")){mode=PER_FILE;}
				else if(b.equalsIgnoreCase("persequence") || b.equalsIgnoreCase("percontig")){mode=PER_SEQUENCE;}
				else if(b.equalsIgnoreCase("perid")){mode=PER_ID;}
				else{throw new RuntimeException("Unknown mode: "+b);}
			}else if(a.equals("perfile")){
				mode=PER_FILE;
			}else if(a.equals("persequence") || a.equals("percontig")){
				mode=PER_SEQUENCE;
			}else if(a.equals("perid")){
				mode=PER_ID;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(Tools.looksLikeInputStream(arg) && inFiles==null){
				inFiles=arg.split(",");
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			Parser.processQuality();
			overwrite=parser.overwrite;
			out=parser.out1;
			if(inFiles==null && parser.in1!=null){
				inFiles=parser.in1.split(",");
			}
		}

		assert(inFiles!=null && inFiles.length>0) : "No input files specified.";
		assert(out!=null) : "No output file specified.";
		if(mode!=PER_FILE){
			throw new RuntimeException("Only perfile mode is currently supported.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();

		for(int f=0; f<inFiles.length; f++){
			final String path=inFiles[f];
			if(verbose){outstream.println("Processing "+path);}

			DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, seed, 0f);

			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			final Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
			cris.start();

			long bases=0, reads=0;
			int contigs=0;
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> list=(ln!=null ? ln.list : null);
			while(ln!=null && list!=null && list.size()>0){
				for(Read r : list){
					ddl.hash(r);
					reads+=r.pairCount();
					bases+=r.pairLength();
					contigs++;
					if(r.mate!=null){contigs++;}
				}
				cris.returnList(ln);
				ln=cris.nextList();
				list=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			ReadWrite.closeStreams(cris);

			DDLRecord rec=new DDLRecord(ddl, f, baseName(path));
			rec.filename=new java.io.File(path).getName();
			rec.bases=bases;
			rec.contigs=contigs;
			rec.cardinality=ddl.cardinality();
			records.add(rec);

			if(verbose){outstream.println("  "+reads+" reads, "+bases+" bases, cardinality="+ddl.cardinality());}
		}

		DDLLoader.writeFile(records, out, overwrite);

		t.stop();
		outstream.println("Wrote "+records.size()+" DDL sketches to "+out);
		outstream.println("Time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Extracts base filename without path or extensions. */
	static String baseName(String path){
		String name=new File(path).getName();
		while(name.contains(".")){
			name=name.substring(0, name.lastIndexOf('.'));
		}
		return name;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String[] inFiles;
	private String out;
	private int k=31;
	private int buckets=2048;
	private long seed=12345L;
	private int mode=PER_FILE;
	private boolean overwrite=false;
	private boolean verbose=false;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int PER_FILE=0, PER_SEQUENCE=1, PER_ID=2;
	private static final PrintStream outstream=System.err;
}
