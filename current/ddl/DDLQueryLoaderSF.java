package ddl;

import java.util.ArrayList;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Single-file DDL query loader with multithreaded hashing and inline gene-calling.
 * Supports two modes:
 *   perFile=false: N threads share one Streamer, each accumulating a partial DDL.
 *     After join, partial DDLs are merged into one DDLRecord.
 *   perFile=true (perContig): N threads share one Streamer, each producing one
 *     DDLRecord per contig.  Contigs shorter than minSketchLength are skipped.
 *
 * Gene-calling behavior:
 *   perFile mode: each thread stops gene-calling after finding its first 16S or 18S.
 *   perContig mode: every contig >= MIN_GENE_CALL_LENGTH is gene-called independently.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
public class DDLQueryLoaderSF {

	/**
	 * Loads a single file in per-file mode (one DDLRecord for the whole file).
	 * @param fname Path to the input sequence file
	 * @param buckets Number of DDL buckets
	 * @param k K-mer length
	 * @param useSSU Whether to do inline gene-calling for SSU
	 * @param maxThreads Maximum worker threads
	 * @return A single DDLRecord representing the entire file
	 */
	public static DDLRecord loadPerFile(String fname, int buckets, int k,
			boolean useSSU, int maxThreads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		final int threads=Math.max(1, maxThreads);
		ProcessThread[] workers=new ProcessThread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new ProcessThread(cris, buckets, k, useSSU, false, 0);
			workers[i].start();
		}
		for(ProcessThread w : workers){
			while(w.getState()!=Thread.State.TERMINATED){
				try{w.join();}catch(InterruptedException e){e.printStackTrace();}
			}
		}
		ReadWrite.closeStreams(cris);

		DynamicDemiLog merged=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
		long bases=0;
		int contigs=0;
		byte[] r16S=null, r18S=null;
		for(ProcessThread w : workers){
			merged.add(w.partialDDL);
			bases+=w.basesT;
			contigs+=w.contigsT;
			if(r16S==null && w.r16S!=null){r16S=w.r16S;}
			if(r18S==null && w.r18S!=null){r18S=w.r18S;}
		}

		DDLRecord rec=new DDLRecord(merged, -1, -1, fname);
		rec.bases=bases;
		rec.contigs=contigs;
		rec.cardinality=merged.cardinality();
		rec.r16S=r16S;
		rec.r18S=r18S;
		return rec;
	}

	/**
	 * Loads a single file in per-contig mode (one DDLRecord per contig).
	 * @param fname Path to the input sequence file
	 * @param buckets Number of DDL buckets
	 * @param k K-mer length
	 * @param useSSU Whether to do inline gene-calling for SSU
	 * @param maxThreads Maximum worker threads
	 * @param minSketchLen Minimum contig length to sketch (shorter contigs skipped)
	 * @return ArrayList of DDLRecords, one per qualifying contig, in input order
	 */
	public static ArrayList<DDLRecord> loadPerContig(String fname, int buckets, int k,
			boolean useSSU, int maxThreads, int minSketchLen){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		final int threads=Math.max(1, maxThreads);
		ProcessThread[] workers=new ProcessThread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new ProcessThread(cris, buckets, k, useSSU, true, minSketchLen);
			workers[i].start();
		}
		for(ProcessThread w : workers){
			while(w.getState()!=Thread.State.TERMINATED){
				try{w.join();}catch(InterruptedException e){e.printStackTrace();}
			}
		}
		ReadWrite.closeStreams(cris);

		ArrayList<long[]> entries=new ArrayList<>();
		for(int t=0; t<workers.length; t++){
			ProcessThread w=workers[t];
			for(int j=0; j<w.perContigRecords.size(); j++){
				entries.add(new long[]{w.numericIDs.get(j), t, j});
			}
		}
		entries.sort((a, b)->Long.compare(a[0], b[0]));

		ArrayList<DDLRecord> results=new ArrayList<>(entries.size());
		for(long[] entry : entries){
			results.add(workers[(int)entry[1]].perContigRecords.get((int)entry[2]));
		}
		return results;
	}

	/**
	 * Gene-calls a file and returns ALL SSU sequences found.
	 * No DDL hashing — just gene-calling.  Each SSU found becomes an SSURecord
	 * with provenance metadata (contig name, start, strand, source file).
	 * @param fname Path to the input sequence file
	 * @param maxThreads Maximum worker threads
	 * @param minGeneCallLen Minimum contig length for gene-calling
	 * @return ArrayList of SSURecords in input order
	 */
	public static ArrayList<SSURecord> findAllSSU(String fname, int maxThreads, int minGeneCallLen){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		final int threads=Math.max(1, maxThreads);
		SSUFindThread[] workers=new SSUFindThread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new SSUFindThread(cris, minGeneCallLen, fname);
			workers[i].start();
		}
		for(SSUFindThread w : workers){
			while(w.getState()!=Thread.State.TERMINATED){
				try{w.join();}catch(InterruptedException e){e.printStackTrace();}
			}
		}
		ReadWrite.closeStreams(cris);

		ArrayList<long[]> entries=new ArrayList<>();
		for(int t=0; t<workers.length; t++){
			SSUFindThread w=workers[t];
			for(int j=0; j<w.ssuRecords.size(); j++){
				entries.add(new long[]{w.ssuRecords.get(j).numericID, t, j});
			}
		}
		entries.sort((a, b)->Long.compare(a[0], b[0]));

		ArrayList<SSURecord> results=new ArrayList<>(entries.size());
		for(long[] entry : entries){
			results.add(workers[(int)entry[1]].ssuRecords.get((int)entry[2]));
		}
		return results;
	}

	/**
	 * Gene-calls pre-parsed Read objects in memory and returns all SSU sequences found.
	 * No file I/O — accepts Read objects directly.
	 * @param reads Pre-parsed Read objects (from HTTP body or other in-memory source)
	 * @param sourceName Label for provenance (e.g. "http_query")
	 * @param minGeneCallLen Minimum contig length for gene-calling
	 * @return ArrayList of SSURecords
	 */
	public static ArrayList<SSURecord> findAllSSU(ArrayList<Read> reads,
			String sourceName, int minGeneCallLen){
		final GeneCaller caller=GeneTools.makeGeneCaller();
		ArrayList<SSURecord> results=new ArrayList<>();
		for(int i=0; i<reads.size(); i++){
			Read r=reads.get(i);
			if(r.bases==null || r.length()<minGeneCallLen){continue;}
			ArrayList<Orf> genes=caller.callGenes(r);
			if(genes==null){continue;}
			for(Orf orf : genes){
				if(!orf.is16S() && !orf.is18S()){continue;}
				SSURecord rec=new SSURecord();
				rec.bases=CallGenes.fetch(orf, r).bases;
				rec.riboType=(orf.is16S() ? DDLRecord.RIBO_16S : DDLRecord.RIBO_18S);
				rec.contigName=r.id;
				rec.start=orf.start;
				rec.strand=(byte)(orf.strand==0 ? '+' : '-');
				rec.numericID=i;
				rec.fileName=sourceName;
				results.add(rec);
			}
		}
		return results;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Thread for finding all SSU sequences via gene-calling.  No DDL hashing. */
	static class SSUFindThread extends Thread {

		SSUFindThread(Streamer cris_, int minGeneCallLen_, String fileName_){
			cris=cris_;
			minGeneCallLen=minGeneCallLen_;
			fileName=fileName_;
		}

		@Override
		public void run(){
			final GeneCaller caller=GeneTools.makeGeneCaller();
			ListNum<Read> ln=cris.nextList();
			while(ln!=null && ln.size()>0){
				for(Read r : ln){
					if(r.bases==null || r.length()<minGeneCallLen){continue;}
					ArrayList<Orf> genes=caller.callGenes(r);
					if(genes==null){continue;}
					for(Orf orf : genes){
						if(!orf.is16S() && !orf.is18S()){continue;}
						SSURecord rec=new SSURecord();
						rec.bases=CallGenes.fetch(orf, r).bases;
						rec.riboType=(orf.is16S() ? DDLRecord.RIBO_16S : DDLRecord.RIBO_18S);
						rec.contigName=r.id;
						rec.start=orf.start;
						rec.strand=(byte)(orf.strand==0 ? '+' : '-');
						rec.numericID=r.numericID;
						rec.fileName=fileName;
						ssuRecords.add(rec);
					}
				}
				cris.returnList(ln);
				ln=cris.nextList();
			}
			cris.returnList(ln);
		}

		final ArrayList<SSURecord> ssuRecords=new ArrayList<>();
		private final Streamer cris;
		private final int minGeneCallLen;
		private final String fileName;
	}

	static class ProcessThread extends Thread {

		ProcessThread(Streamer cris_, int buckets_, int k_, boolean useSSU_,
				boolean perContig_, int minSketchLen_){
			cris=cris_;
			buckets=buckets_;
			k=k_;
			useSSU=useSSU_;
			perContig=perContig_;
			minSketchLen=minSketchLen_;
		}

		@Override
		public void run(){
			final GeneCaller caller=(useSSU ? GeneTools.makeGeneCaller() : null);
			boolean needSSU=useSSU;

			if(!perContig){
				partialDDL=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
			}

			ListNum<Read> ln=cris.nextList();
			while(ln!=null && ln.size()>0){
				for(Read r : ln){
					if(r.bases==null){continue;}
					final int len=r.length();

					if(!perContig){
						partialDDL.hash(r);
						basesT+=len;
						contigsT++;
						if(r.mate!=null){basesT+=r.mate.length();}
						if(needSSU && len>=MIN_GENE_CALL_LENGTH){
							ArrayList<Orf> genes=caller.callGenes(r);
							if(genes!=null){
								for(Orf orf : genes){
									if(orf.is16S()){r16S=CallGenes.fetch(orf, r).bases; needSSU=false; break;}
									else if(orf.is18S()){r18S=CallGenes.fetch(orf, r).bases; needSSU=false; break;}
								}
							}
						}
					}else{
						if(len<minSketchLen){continue;}
						DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
						ddl.hash(r);
						DDLRecord rec=new DDLRecord(ddl, -1, -1, r.id);
						rec.bases=len;
						rec.contigs=1;
						rec.cardinality=ddl.cardinality();
						if(caller!=null && len>=MIN_GENE_CALL_LENGTH){
							ArrayList<Orf> genes=caller.callGenes(r);
							if(genes!=null){
								for(Orf orf : genes){
									if(orf.is16S()){rec.r16S=CallGenes.fetch(orf, r).bases; break;}
									else if(orf.is18S()){rec.r18S=CallGenes.fetch(orf, r).bases; break;}
								}
							}
						}
						perContigRecords.add(rec);
						numericIDs.add(r.numericID);
					}
				}
				cris.returnList(ln);
				ln=cris.nextList();
			}
			cris.returnList(ln);
		}

		/* Per-file mode fields */
		DynamicDemiLog partialDDL;
		long basesT;
		int contigsT;
		byte[] r16S;
		byte[] r18S;

		/* Per-contig mode fields */
		final ArrayList<DDLRecord> perContigRecords=new ArrayList<>();
		final ArrayList<Long> numericIDs=new ArrayList<>();

		/* Config */
		private final Streamer cris;
		private final int buckets, k;
		private final boolean useSSU, perContig;
		private final int minSketchLen;
	}

	static final int MIN_GENE_CALL_LENGTH=800;
}
