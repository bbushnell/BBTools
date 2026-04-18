package ddl;

import java.io.File;

import cardinality.DynamicDemiLog;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import HashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
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
 * as A48-encoded TSV. Supports perfile (with multithreading), persequence,
 * and pertid modes.
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
				else if(b.equalsIgnoreCase("pertid") || b.equalsIgnoreCase("pertaxid")){mode=PER_TID;}
				else{throw new RuntimeException("Unknown mode: "+b);}
			}else if(a.equals("perfile")){
				mode=PER_FILE;
			}else if(a.equals("persequence") || a.equals("percontig")){
				mode=PER_SEQUENCE;
			}else if(a.equals("pertid") || a.equals("pertaxid")){
				mode=PER_TID;
			}else if(a.equals("parsetaxid") || a.equals("parsetid")){
				parseTaxid=Parse.parseBoolean(b);
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
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		if(mode==PER_FILE){processPerFile(t);}
		else if(mode==PER_SEQUENCE){processPerSequence(t);}
		else if(mode==PER_TID){processPerTid(t);}
		else{throw new RuntimeException("Unknown mode: "+mode);}
	}

	void processPerFile(Timer t){
		final int threads=Tools.min(Shared.threads(), inFiles.length);
		if(threads>1){processPerFileMulti(t, threads); return;}

		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();
		for(int f=0; f<inFiles.length; f++){
			DDLRecord rec=processOneFile(f);
			if(rec!=null){records.add(rec);}
		}
		DDLLoader.writeFile(records, out, overwrite, k, seed);

		t.stop();
		outstream.println("Wrote "+records.size()+" DDL sketches to "+out);
		outstream.println("Time: \t"+t);
	}

	void processPerFileMulti(Timer t, int threads){
		final ArrayList<DDLRecord> allRecords=
			Collections.synchronizedList(new ArrayList<DDLRecord>());
		final AtomicInteger fileIndex=new AtomicInteger(0);

		Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new Thread("DDLWriter-"+i){
				@Override
				public void run(){
					for(int f=fileIndex.getAndIncrement(); f<inFiles.length;
						f=fileIndex.getAndIncrement()){
						DDLRecord rec=processOneFile(f);
						if(rec!=null){allRecords.add(rec);}
					}
				}
			};
			workers[i].start();
		}
		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}

		//Sort by id (load order) for deterministic output
		ArrayList<DDLRecord> records=new ArrayList<>(allRecords);
		records.sort((a, b) -> a.id-b.id);
		DDLLoader.writeFile(records, out, overwrite, k, seed);

		t.stop();
		outstream.println("Wrote "+records.size()+" DDL sketches to "+out+" using "+threads+" threads");
		outstream.println("Time: \t"+t);
	}

	/** Process a single file, returning one DDLRecord. Thread-safe. */
	private DDLRecord processOneFile(int f){
		final String path=inFiles[f];
		if(verbose){outstream.println("Processing "+path);}

		DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, seed, 0f);

		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		final Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		long bases=0;
		long gcBases=0, atBases=0;
		int contigs=0;
		int taxID=-1;
		String firstName=null;
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> list=(ln!=null ? ln.list : null);
		while(ln!=null && list!=null && list.size()>0){
			for(Read r : list){
				ddl.hash(r);
				bases+=r.pairLength();
				contigs++;
				if(r.mate!=null){contigs++;}
				if(firstName==null){firstName=r.id;}
				if(r.bases!=null){
					for(byte b : r.bases){
						int x=AminoAcid.baseToNumber[b];
						if(x==1 || x==2){gcBases++;}
						else if(x==0 || x==3){atBases++;}
					}
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			list=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);

		if(parseTaxid){
			taxID=parseTID(new File(path).getName());
			if(taxID<1 && firstName!=null){taxID=parseTID(firstName);}
		}

		DDLRecord rec=new DDLRecord(ddl, f, taxID, baseName(path));
		rec.filename=new File(path).getName();
		rec.bases=bases;
		rec.contigs=contigs;
		rec.cardinality=ddl.cardinality();
		if(gcBases+atBases>0){rec.gc=gcBases*1f/(gcBases+atBases);}
		return rec;
	}

	void processPerSequence(Timer t){
		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();
		int nextId=0;

		for(int f=0; f<inFiles.length; f++){
			final String path=inFiles[f];
			if(verbose){outstream.println("Processing "+path);}

			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			final Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
			cris.start();

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> list=(ln!=null ? ln.list : null);
			while(ln!=null && list!=null && list.size()>0){
				for(Read r : list){
					DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, seed, 0f);
					ddl.hash(r);

					int taxID=-1;
					if(parseTaxid){
						taxID=parseTID(r.id);
						assert(taxID>0) : "parseTaxid=true but no TID in header: "+r.id;
					}

					DDLRecord rec=new DDLRecord(ddl, nextId++, taxID, r.id);
					rec.bases=r.length();
					rec.contigs=1;
					rec.cardinality=ddl.cardinality();
					if(r.bases!=null){
						long gc=0, at=0;
						for(byte b : r.bases){
							int x=AminoAcid.baseToNumber[b];
							if(x==1 || x==2){gc++;}
							else if(x==0 || x==3){at++;}
						}
						if(gc+at>0){rec.gc=gc*1f/(gc+at);}
					}
					records.add(rec);
				}
				cris.returnList(ln);
				ln=cris.nextList();
				list=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			ReadWrite.closeStreams(cris);
		}

		DDLLoader.writeFile(records, out, overwrite, k, seed);

		t.stop();
		outstream.println("Wrote "+records.size()+" DDL sketches to "+out);
		outstream.println("Time: \t"+t);
	}

	void processPerTid(Timer t){
		final HashMap<Integer, DDLRecord> tidMap=new HashMap<>();
		final HashMap<Integer, long[]> gcMap=new HashMap<>(); //{gcBases, atBases}
		int nextId=0;

		for(int f=0; f<inFiles.length; f++){
			final String path=inFiles[f];
			if(verbose){outstream.println("Processing "+path);}

			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			final Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
			cris.start();

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> list=(ln!=null ? ln.list : null);
			while(ln!=null && list!=null && list.size()>0){
				for(Read r : list){
					int tid=parseTID(r.id);
					if(tid<1){
						outstream.println("WARNING: No TID in header, skipping: "+r.id);
						continue;
					}
					DDLRecord rec=tidMap.get(tid);
					if(rec==null){
						DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, seed, 0f);
						rec=new DDLRecord(ddl, nextId++, tid, r.id);
						tidMap.put(tid, rec);
						gcMap.put(tid, new long[2]);
					}
					rec.ddl.hash(r);
					rec.bases+=r.length();
					rec.contigs++;
					if(r.bases!=null){
						long[] gc=gcMap.get(tid);
						for(byte b : r.bases){
							int x=AminoAcid.baseToNumber[b];
							if(x==1 || x==2){gc[0]++;}
							else if(x==0 || x==3){gc[1]++;}
						}
					}
				}
				cris.returnList(ln);
				ln=cris.nextList();
				list=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			ReadWrite.closeStreams(cris);
		}

		//Finalize and sort
		ArrayList<DDLRecord> records=new ArrayList<>(tidMap.values());
		for(DDLRecord rec : records){
			rec.cardinality=rec.ddl.cardinality();
			long[] gc=gcMap.get(rec.taxID);
			if(gc!=null && gc[0]+gc[1]>0){rec.gc=gc[0]*1f/(gc[0]+gc[1]);}
		}
		Collections.sort(records);
		DDLLoader.writeFile(records, out, overwrite, k, seed);

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

	/** Parse TID from header or filename like "tid|504728|..." or "tid_504728_..." */
	static int parseTID(String s){
		if(s==null){return -1;}
		int pos=s.indexOf("tid|");
		if(pos<0){pos=s.indexOf("tid_");}
		if(pos<0){return -1;}
		long id=0;
		for(int i=pos+4; i<s.length(); i++){
			char c=s.charAt(i);
			if(c<'0' || c>'9'){break;}
			id=id*10+(c-'0');
		}
		return (id>0 && id<Integer.MAX_VALUE) ? (int)id : -1;
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
	private boolean parseTaxid=true;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int PER_FILE=0, PER_SEQUENCE=1, PER_TID=2;
	private static final PrintStream outstream=System.err;
}
