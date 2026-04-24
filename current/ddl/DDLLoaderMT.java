package ddl;

import java.util.ArrayList;

import cardinality.DynamicDemiLog;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;
import stream.JobQueue;
import structures.ListNum;

/**
 * Multithreaded DDL file loader using producer-consumer pattern.
 * Single producer reads lines and bundles complete DDL records;
 * N worker threads parse records into DDLRecord objects in parallel.
 *
 * @author Eru
 * @date April 2026
 */
public class DDLLoaderMT {

	public static ArrayList<DDLRecord> loadFile(String path, int k){
		return loadFile(path, k, Shared.threads());
	}

	public static ArrayList<DDLRecord> loadFile(String path, int k, int threads){
		if(threads<2){return DDLLoader.loadFile(path, k);}

		final FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);
		final ByteFile bf=ByteFile.makeByteFile(ff);

		final JobQueue<ListNum<ArrayList<byte[]>>> queue=
				new JobQueue<>(threads+4, true, true, 0);

		@SuppressWarnings("unchecked")
		final ArrayList<DDLRecord>[] allResults=new ArrayList[threads];
		final Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			final ArrayList<DDLRecord> threadResults=new ArrayList<>();
			allResults[i]=threadResults;
			final int kk=k;
			workers[i]=new Thread("DDLLoaderMT-"+i){
				@Override
				public void run(){
					consume(queue, threadResults, kk);
				}
			};
			workers[i].start();
		}

		produce(bf, queue);

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}
		bf.close();

		final ArrayList<DDLRecord> merged=new ArrayList<>();
		for(ArrayList<DDLRecord> results : allResults){
			merged.addAll(results);
		}
		return merged;
	}

	private static void produce(ByteFile bf,
			JobQueue<ListNum<ArrayList<byte[]>>> queue){
		long id=0;
		ArrayList<ArrayList<byte[]>> bundle=new ArrayList<>(RECORDS_PER_BUNDLE);
		ArrayList<byte[]> currentRecord=new ArrayList<>(12);

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if((Tools.startsWith(line, "#tid") || Tools.startsWith(line, "#id"))
					&& !currentRecord.isEmpty()){
				bundle.add(currentRecord);
				currentRecord=new ArrayList<>(12);
				currentRecord.add(line);
				if(bundle.size()>=RECORDS_PER_BUNDLE){
					queue.add(new ListNum<>(bundle, id++));
					bundle=new ArrayList<>(RECORDS_PER_BUNDLE);
				}
			}else{
				currentRecord.add(line);
			}
		}

		if(!currentRecord.isEmpty()){bundle.add(currentRecord);}
		if(!bundle.isEmpty()){
			queue.add(new ListNum<>(bundle, id++, false, true));
		}else{
			queue.poison(new ListNum<>(null, id, true, false), true);
		}
	}

	private static void consume(JobQueue<ListNum<ArrayList<byte[]>>> queue,
			ArrayList<DDLRecord> results, int k){
		final LineParser1 lp=new LineParser1('\t');
		for(ListNum<ArrayList<byte[]>> batch=queue.take();
				batch!=null; batch=queue.take()){
			for(ArrayList<byte[]> recordLines : batch){
				DDLRecord rec=parseRecord(recordLines, lp, k);
				if(rec!=null){results.add(rec);}
			}
		}
	}

	private static DDLRecord parseRecord(ArrayList<byte[]> lines, LineParser1 lp, int k){
		long recId=-1;
		int tid=-1;
		String name=null, file=null;
		long bases=0;
		int contigs=0;
		float gc=-1;
		byte[] dataLine=null;

		for(byte[] line : lines){
			if(line.length<1){continue;}
			if(line[0]=='#'){
				lp.set(line);
				if(lp.terms()<2){continue;}
				if(lp.termEquals("#id", 0)){recId=lp.parseLong(1);}
				else if(lp.termEquals("#tid", 0)){tid=(int)lp.parseLong(1);}
				else if(lp.termEquals("#name", 0)){name=lp.parseString(1);}
				else if(lp.termEquals("#file", 0)){file=lp.parseString(1);}
				else if(lp.termEquals("#bases", 0)){bases=lp.parseLong(1);}
				else if(lp.termEquals("#contigs", 0)){contigs=(int)lp.parseLong(1);}
				else if(lp.termEquals("#gc", 0)){gc=lp.parseFloat(1);}
			}else{
				dataLine=line;
			}
		}

		if(dataLine==null){return null;}

		DynamicDemiLog ddl=DDLLoader.parseDDL(dataLine, lp, k);
		DDLRecord rec=new DDLRecord(ddl, recId, tid, name);
		rec.filename=file;
		rec.bases=bases;
		rec.contigs=contigs;
		rec.gc=gc;
		rec.cardinality=ddl.cardinality();
		return rec;
	}

	private static final int RECORDS_PER_BUNDLE=8;

}
