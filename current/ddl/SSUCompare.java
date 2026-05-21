package ddl;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import parse.Parser;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import prok.ProkObject;
import server.ServerTools;
import shared.Resources;
import shared.Shared;
import structures.StringNum;
import shared.Timer;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;

/**
 * SSU-specific comparison tool.  Compares 16S/18S ribosomal sequences
 * against a pre-built SSU reference database using DDL sketching.
 *
 * Two modes:
 *   Default: input sequences are assumed to be SSUs.  Each is classified
 *     as 16S or 18S by alignment against consensus sequences.
 *   Call mode: input is genomic sequence.  Gene-calling finds SSUs,
 *     potentially multiple per contig.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
public class SSUCompare {

	public static void main(String[] args){
		if(args.length<1){
			System.err.println("Usage: SSUCompare <ssu1.fa> [ssu2.fa ...] ref=<ssu_ddls.tsv> [records=5]");
			System.err.println("   or: SSUCompare <genome.fa> call ref=<ssu_ddls.tsv>");
			System.exit(1);
		}

		int k=19, buckets=128, maxRecords=5, minHits=8, buffer=0;
		boolean useIndex=true, callMode=false, alignSSU=true, banSelf=false;
		boolean local=false, loud=false;
		String address=null;
		String refFile=null, ref16sFile=null, ref18sFile=null, queryFile=null;
		ArrayList<String> inFiles=new ArrayList<>();
		DDLFormatter formatter=new DDLFormatter();
		Parser parser=new Parser();
		DynamicDemiLog.setExponent(4);

		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(formatter.parse(args[i], a, b)){/* handled */}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("ref16s") || a.equals("ref16")){ref16sFile=b;}
			else if(a.equals("ref18s") || a.equals("ref18")){ref18sFile=b;}
			else if(a.equals("qf") || a.equals("queryfile")){queryFile=b;}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("buffer")){buffer=Integer.parseInt(b);}
			else if(a.equals("index")){useIndex=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("call") || a.equals("callssu")){callMode=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("align") || a.equals("alignssu")){alignSSU=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("banself")){banSelf=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("loud")){loud=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("address")){address=b;}
			else if(a.equals("local")){local=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("in")){inFiles.add(b);}
			else if(parser.parse(args[i], a, b)){/* handled */}
			else if(b==null && !a.startsWith("-")){inFiles.add(args[i]);}
		}

		int threads=Shared.threads();
		if(!local && address==null && refFile==null && ref16sFile==null && ref18sFile==null){
			address=DEFAULT_ADDRESS;
		}
		if(address!=null && !local){
			sendToServer(inFiles, address, callMode, maxRecords, minHits, buffer, formatter);
			return;
		}

		if(refFile==null && ref16sFile==null && ref18sFile==null){
			refFile=Resources.find("?ssuSketchDDL.tsv.gz");
		}
		if(queryFile==null && inFiles.isEmpty()){
			System.err.println("ERROR: no input files or queryfile specified.");
			System.exit(1);
		}

		formatter.useAlignmentANI=true;
		formatter.printSSU=false;
		formatter.printBases=false;
		formatter.printCompleteness=false;
		formatter.printQueryName=false;
		formatter.printRank=true;
		formatter.printType=true;
		formatter.printQLen=true;
		formatter.printRLen=true;
		formatter.printFile=true;
		formatter.printContig=true;
		formatter.printStart=true;
		formatter.printStrand=true;

		run(inFiles, queryFile, refFile, ref16sFile, ref18sFile, k, buckets, maxRecords, buffer, minHits, useIndex, callMode, alignSSU, banSelf, threads, loud, formatter);
	}

	private static void run(ArrayList<String> inFiles, String queryFile, String refPath, String ref16sPath, String ref18sPath,
			int k, int buckets, int maxRecords, int buffer, int minHits, boolean useIndex, boolean callMode,
			boolean alignSSU, boolean banSelf, int threads, boolean loud, DDLFormatter formatter){
		Timer t=new Timer();
		if(loud){System.err.println("SSUCompare: "+(callMode ? "call" : "ssu")+" mode  threads="+threads
			+"  k="+k+"  exponent="+DynamicDemiLog.exponentBits()+"  buckets="+buckets);}

		/* --- Load refs --- */

		long ts=System.nanoTime();
		ArrayList<DDLRecord> refs=new ArrayList<>();
		map.IntObjectMap<byte[]>[] maps=DDLSSULoader.loadSSUMapsDefaults();
		map.IntObjectMap<byte[]> map16=maps[0], map18=maps[1];

		if(ref16sPath!=null || ref18sPath!=null){
			if(ref16sPath!=null){
				System.err.println("Loading 16S references from "+ref16sPath+"...");
				ArrayList<DDLRecord> r16=DDLLoaderMT.loadFile(ref16sPath, k, threads);
				for(DDLRecord rec : r16){
					if(rec.taxID>0 && map16!=null){
						byte[] seq=map16.get(rec.taxID);
						if(seq!=null){rec.r16S=seq;}
					}
				}
				System.err.println("Loaded "+r16.size()+" 16S refs.");
				refs.addAll(r16);
			}
			if(ref18sPath!=null){
				System.err.println("Loading 18S references from "+ref18sPath+"...");
				ArrayList<DDLRecord> r18=DDLLoaderMT.loadFile(ref18sPath, k, threads);
				for(DDLRecord rec : r18){
					if(rec.taxID>0 && map18!=null){
						byte[] seq=map18.get(rec.taxID);
						if(seq!=null){rec.r18S=seq;}
					}
				}
				System.err.println("Loaded "+r18.size()+" 18S refs.");
				refs.addAll(r18);
			}
		}else{
			System.err.println("Loading SSU references from "+refPath+"...");
			refs=DDLLoaderMT.loadFile(refPath, k, threads);
			int a16=0, a18=0;
			for(DDLRecord rec : refs){
				if(rec.taxID<=0){continue;}
				boolean is16S=(rec.filename!=null && rec.filename.contains("16S"));
				boolean is18S=(rec.filename!=null && rec.filename.contains("18S"));
				if(is16S && map16!=null){
					byte[] seq=map16.get(rec.taxID);
					if(seq!=null){rec.r16S=seq; a16++;}
				}else if(is18S && map18!=null){
					byte[] seq=map18.get(rec.taxID);
					if(seq!=null){rec.r18S=seq; a18++;}
				}else if(map16!=null){
					byte[] seq=map16.get(rec.taxID);
					if(seq!=null){rec.r16S=seq; a16++;}
					else if(map18!=null){
						seq=map18.get(rec.taxID);
						if(seq!=null){rec.r18S=seq; a18++;}
					}
				}
			}
			System.err.println("Attached "+a16+" 16S and "+a18+" 18S sequences to "+refs.size()+" refs.");
		}
		long tRefLoad=System.nanoTime()-ts;
		System.err.println("Total: "+refs.size()+" reference SSU DDLs in "+fmt(tRefLoad)+" seconds.");
		long tSSULoad=0;

		DDLIndex index=null;
		if(useIndex){
			ts=System.nanoTime();
			index=new DDLIndex(refs.get(0).ddl.buckets);
			index.addAll(refs, threads);
			long tIndex=System.nanoTime()-ts;
			System.err.println("Built index in "+fmt(tIndex)+" seconds.");
		}

		/* --- Load queries --- */

		ts=System.nanoTime();
		ArrayList<DDLRecord> queries;
		if(queryFile!=null){
			System.err.println("Loading pre-built queries from "+queryFile+"...");
			queries=DDLLoaderMT.loadFile(queryFile, k, threads);
			if(alignSSU){
				int qa16=0, qa18=0;
				for(DDLRecord rec : queries){
					if(rec.taxID<=0){continue;}
					boolean is16S=(rec.filename!=null && rec.filename.contains("16S"));
					boolean is18S=(rec.filename!=null && rec.filename.contains("18S"));
					if(is16S && map16!=null){
						byte[] seq=map16.get(rec.taxID);
						if(seq!=null){rec.r16S=seq; qa16++;}
					}else if(is18S && map18!=null){
						byte[] seq=map18.get(rec.taxID);
						if(seq!=null){rec.r18S=seq; qa18++;}
					}else if(map16!=null){
						byte[] seq=map16.get(rec.taxID);
						if(seq!=null){rec.r16S=seq; qa16++;}
						else if(map18!=null){
							seq=map18.get(rec.taxID);
							if(seq!=null){rec.r18S=seq; qa18++;}
						}
					}
				}
				System.err.println("Attached "+qa16+" 16S and "+qa18+" 18S to "+queries.size()+" queries.");
			}
		}else if(callMode){
			queries=loadCallMode(inFiles, k, buckets, threads);
		}else{
			queries=loadSSUMode(inFiles, k, buckets);
		}
		long tQueryLoad=System.nanoTime()-ts;
		System.err.println("Loaded "+queries.size()+" SSU queries in "+fmt(tQueryLoad)+" seconds.");

		if(queries.isEmpty()){
			System.err.println("No SSU queries found.");
			t.stop();
			System.err.println("Total time: \t"+t);
			return;
		}

		/* --- Compare --- */

		final int nRefs=refs.size();
		final int nQueries=queries.size();

		@SuppressWarnings("unchecked")
		final ArrayList<DDLComparison>[] allResults=new ArrayList[nQueries];
		final AtomicLong nextQuery=new AtomicLong(0);
		final AtomicLong totalComparisons=new AtomicLong(0);

		ts=System.nanoTime();
		CompareThread[] workers=new CompareThread[threads];
		for(int wi=0; wi<threads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery, allResults,
				nQueries, nRefs, k, maxRecords, buffer, minHits, useIndex, index,
				totalComparisons, alignSSU, banSelf);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}
		long tCompare=System.nanoTime()-ts;

		long totalAlignments=0;
		for(CompareThread w : workers){totalAlignments+=w.alignCountT;}
		System.err.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+"s"
			+"  ("+totalAlignments+" SSU alignments inline)");

		/* --- Format --- */

		ts=System.nanoTime();
		ByteBuilder bb=new ByteBuilder();
		boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
		if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
		for(int qi=0; qi<nQueries; qi++){
			if(allResults[qi]==null){continue;}
			int rank=0;
			for(DDLComparison c : allResults[qi]){
				if(c.equal<minHits){continue;}
				c.rank=++rank;
				formatter.format(c, bb);
			}
		}
		if(json){formatter.jsonEnd(bb);}
		System.out.print(bb);
		long tFormat=System.nanoTime()-ts;

		t.stop();
		if(loud){
			System.err.println("\nSubphase timing:");
			System.err.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
			System.err.println("  SSU attach:    \t"+fmt(tSSULoad)+"s");
			System.err.println("  Query load:    \t"+fmt(tQueryLoad)+"s  ("+nQueries+" queries)");
			System.err.println("  Compare+align: \t"+fmt(tCompare)+"s  ("+totalAlignments+" alignments)");
			System.err.println("  Format:        \t"+fmt(tFormat)+"s");
		}
		System.err.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Mode 1: SSU Input        ----------------*/
	/*--------------------------------------------------------------*/

	private static ArrayList<DDLRecord> loadSSUMode(ArrayList<String> inFiles, int k, int buckets){
		byte[] consensus16S=loadFirstSequence("?16S_consensus_sequence.fa");
		byte[] consensus18S=loadFirstSequence("?18S_consensus_sequence.fa");
		QuantumAligner aligner=new QuantumAligner();
		ArrayList<DDLRecord> queries=new ArrayList<>();

		for(String fname : inFiles){
			FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			ArrayList<Read> reads=StreamerFactory.getReads(-1, false, ff, null, null, null);
			for(Read r : reads){
				if(r.bases==null || r.length()<50){continue;}
				float id16=(consensus16S!=null ? aligner.align(r.bases, consensus16S) : 0);
				float id18=(consensus18S!=null ? aligner.align(r.bases, consensus18S) : 0);
				boolean is16S=(id16>=id18);
				if(id16<0.3f && id18<0.3f){
					System.err.println("Warning: low consensus identity for "+r.id
						+" (16S="+String.format("%.3f", id16)+", 18S="+String.format("%.3f", id18)+")");
				}

				DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
				ddl.hash(r);
				DDLRecord rec=new DDLRecord(ddl, -1, -1, r.id);
				rec.bases=r.length();
				rec.contigs=1;
				rec.cardinality=ddl.cardinality();
				rec.filename=new File(fname).getName();
				rec.contigName=r.id;
				rec.ssuStart=0;
				rec.ssuStrand=(byte)'+';
				if(is16S){rec.r16S=r.bases;}else{rec.r18S=r.bases;}
				queries.add(rec);
			}
		}
		return queries;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Mode 2: Call Mode         ----------------*/
	/*--------------------------------------------------------------*/

	private static ArrayList<DDLRecord> loadCallMode(ArrayList<String> inFiles, int k, int buckets, int threads){
		ProkObject.callCDS=false;
		ProkObject.calltRNA=false;
		ProkObject.call23S=false;
		ProkObject.call5S=false;
		GeneTools.loadPGM();

		ArrayList<DDLRecord> queries=new ArrayList<>();
		for(String fname : inFiles){
			ArrayList<SSURecord> ssus=DDLQueryLoaderSF.findAllSSU(fname, threads, MIN_GENE_CALL_LENGTH);
			System.err.println("Found "+ssus.size()+" SSU sequences in "+fname);
			for(SSURecord ssu : ssus){
				String cname=ssu.contigName;
				if(cname!=null && cname.indexOf(' ')>0){cname=cname.substring(0, cname.indexOf(' '));}
				Read r=new Read(ssu.bases, null, cname+":"+ssu.start+(char)ssu.strand, 0);
				DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
				ddl.hash(r);
				DDLRecord rec=new DDLRecord(ddl, -1, -1, r.id);
				rec.bases=ssu.bases.length;
				rec.contigs=1;
				rec.cardinality=ddl.cardinality();
				rec.filename=new File(ssu.fileName).getName();
				rec.contigName=cname;
				rec.ssuStart=ssu.start;
				rec.ssuStrand=ssu.strand;
				if(ssu.is16S()){rec.r16S=ssu.bases;}else{rec.r18S=ssu.bases;}
				queries.add(rec);
			}
		}
		return queries;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	/** Attaches SSU bytes to refs, ensuring each ref gets only ONE type.
	 * A ref from the 16S file gets r16S only; a ref from 18S gets r18S only.
	 * If a TaxID appears in both files, each ref entry retains only its type. */
	private static void attachSSUSeparate(ArrayList<DDLRecord> refs){
		map.IntObjectMap<byte[]>[] maps=DDLSSULoader.loadSSUMapsDefaults();
		map.IntObjectMap<byte[]> map16=maps[0], map18=maps[1];
		int a16=0, a18=0;
		for(DDLRecord rec : refs){
			if(rec.taxID<=0){continue;}
			byte[] s16=(map16!=null ? map16.get(rec.taxID) : null);
			byte[] s18=(map18!=null ? map18.get(rec.taxID) : null);
			if(s16!=null && s18!=null){
				if(rec.r16S==null && rec.r18S==null){
					rec.r16S=s16; a16++;
				}
			}else if(s16!=null){
				rec.r16S=s16; a16++;
			}else if(s18!=null){
				rec.r18S=s18; a18++;
			}
		}
		System.err.println("Attached "+a16+" 16S and "+a18+" 18S sequences to "+refs.size()+" refs.");
	}

	private static byte[] loadFirstSequence(String path){
		String resolved=dna.Data.findPath(path);
		if(resolved==null){System.err.println("Warning: consensus file not found: "+path); return null;}
		FileFormat ff=FileFormat.testInput(resolved, FileFormat.FASTA, null, true, true);
		ArrayList<Read> reads=StreamerFactory.getReads(-1, false, ff, null, null, null);
		return (reads!=null && !reads.isEmpty()) ? reads.get(0).bases : null;
	}

	private static String fmt(long nanos){return String.format("%.3f", nanos*1e-9);}

	/*--------------------------------------------------------------*/
	/*----------------       Client Mode            ----------------*/
	/*--------------------------------------------------------------*/

	private static void sendToServer(ArrayList<String> inFiles, String address, boolean callMode,
			int maxRecords, int minHits, int buffer, DDLFormatter formatter){
		ByteBuilder bb=new ByteBuilder();
		bb.append("//records=").append(maxRecords).append('\n');
		bb.append("//minhits=").append(minHits).append('\n');
		if(buffer>0){bb.append("//buffer=").append(buffer).append('\n');}
		if(formatter.printLineage){bb.append("//lineage=t\n");}
		if(formatter.printRank){bb.append("//rank=t\n");}
		if(formatter.format==DDLFormatter.FORMAT_JSON){bb.append("//json=t\n");}
		if(callMode){
			ProkObject.callCDS=false;
			ProkObject.calltRNA=false;
			ProkObject.call23S=false;
			ProkObject.call5S=false;
			GeneTools.loadPGM();
			int totalSSU=0;
			for(String fname : inFiles){
				ArrayList<SSURecord> ssus=DDLQueryLoaderSF.findAllSSU(fname, 1, MIN_GENE_CALL_LENGTH);
				System.err.println("Found "+ssus.size()+" SSU sequences in "+fname);
				for(SSURecord ssu : ssus){
					String cname=ssu.contigName;
					if(cname!=null && cname.indexOf(' ')>0){cname=cname.substring(0, cname.indexOf(' '));}
					bb.append('>').append(cname).append(':').append(ssu.start).append((char)ssu.strand).append('\n');
					bb.append(ssu.bases).append('\n');
				}
				totalSSU+=ssus.size();
			}
			System.err.println("Gene-called "+totalSSU+" SSUs locally, sending to server.");
		}else{
			for(String fname : inFiles){
				try{
					byte[] raw=ReadWrite.readRaw(fname);
					if(raw!=null){bb.append(raw);}
				}catch(java.io.IOException e){
					System.err.println("ERROR reading "+fname+": "+e.getMessage());
					System.exit(1);
				}
			}
		}
		if(bb.length()==0){
			System.err.println("ERROR: no input data to send.");
			System.exit(1);
		}
		System.err.println("Sending "+bb.length()+" bytes to "+address);
		StringNum result=ServerTools.sendAndReceive(bb.toBytes(), address);
		if(result!=null && result.s!=null){
			System.out.print(result.s);
		}else{
			System.err.println("ERROR: No response from server at "+address);
			System.exit(1);
		}
	}

	private static final int MIN_GENE_CALL_LENGTH=800;
	private static final String DEFAULT_ADDRESS="https://bbmapservers.jgi.doe.gov/sendclade/findssu/";
}
