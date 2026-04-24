package clade;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import bin.AdjustEntropy;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.JobQueue;
import structures.ListNum;

/**
 * Multithreaded clade spectra loader using producer-consumer pattern.
 * Single producer reads lines and bundles complete records;
 * N worker threads parse records into Clade objects in parallel.
 *
 * @author Eru
 * @date April 2026
 */
public class CladeLoaderMT {

	/*--------------------------------------------------------------*/
	/*----------------        Loading Methods       ----------------*/
	/*--------------------------------------------------------------*/

	static void loadClades(ByteFile bf, ConcurrentHashMap<Integer, Clade> map){
		loadClades(bf, map, loadThreads>0 ? loadThreads : Shared.threads());
	}

	static void loadClades(ByteFile bf, ConcurrentHashMap<Integer, Clade> map, int threads){
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
		if(threads<2){
			CladeLoader.loadClades(bf, map);
			return;
		}

		final JobQueue<ListNum<ArrayList<byte[]>>> queue=
				new JobQueue<>(threads+4, true, true, 0);

		@SuppressWarnings("unchecked")
		final ArrayList<Clade>[] allResults=new ArrayList[threads];
		final Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			final ArrayList<Clade> threadResults=new ArrayList<>();
			allResults[i]=threadResults;
			workers[i]=new Thread("CladeLoaderMT-"+i){
				@Override
				public void run(){
					consume(queue, threadResults);
				}
			};
			workers[i].start();
		}

		produce(bf, queue);

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}

		for(ArrayList<Clade> results : allResults){
			for(Clade c : results){
				Integer key=c.taxID;
				Clade old=map.get(key);
				if(old==null){
					map.put(c.taxID, c);
				}else if(c.bases>0){
					if(CladeLoader.mergeDuplicateTaxIDs){
						old.add(c);
					}else if(c.bases>old.bases){
						map.put(c.taxID, c);
					}
				}
			}
		}
	}

	private static void produce(ByteFile bf,
			JobQueue<ListNum<ArrayList<byte[]>>> queue){
		long id=0;
		ArrayList<ArrayList<byte[]>> bundle=new ArrayList<>(RECORDS_PER_BUNDLE);
		ArrayList<byte[]> currentRecord=new ArrayList<>(20);

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(Tools.startsWith(line, '#') && currentRecord.size()>5){
				bundle.add(currentRecord);
				currentRecord=new ArrayList<>(20);
				currentRecord.add(line);
				if(bundle.size()>=RECORDS_PER_BUNDLE){
					queue.add(new ListNum<>(bundle, id++));
					bundle=new ArrayList<>(RECORDS_PER_BUNDLE);
				}
			}else{
				currentRecord.add(line);
			}
		}

		if(currentRecord.size()>5){bundle.add(currentRecord);}
		if(!bundle.isEmpty()){
			queue.add(new ListNum<>(bundle, id++, false, true));
		}else{
			queue.poison(new ListNum<>(null, id, true, false), true);
		}
	}

	private static void consume(JobQueue<ListNum<ArrayList<byte[]>>> queue,
			ArrayList<Clade> results){
		final LineParser1 lp=new LineParser1('\t');
		for(ListNum<ArrayList<byte[]>> batch=queue.take();
				batch!=null; batch=queue.take()){
			for(ArrayList<byte[]> record : batch){
				if(record.size()>5){
					Clade c=Clade.parseClade(record, lp);
					if(c!=null){results.add(c);}
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Benchmark Main        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		String path=null;
		int threads=Shared.threads();
		int mode=2;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("ref")){path=b;}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("mode")){mode=Integer.parseInt(b);}
		}
		if(path==null){
			System.err.println("Usage: CladeLoaderMT in=<spectra.gz> threads=N mode=0|1|2");
			System.err.println("  mode=0: read lines only (I/O baseline)");
			System.err.println("  mode=1: produce+queue, workers discard (queue overhead)");
			System.err.println("  mode=2: full parse (production path)");
			return;
		}

		Shared.setThreads(threads);
		FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);

		Timer t=new Timer(System.err, false);
		if(mode==0){
			benchReadOnly(ff);
		}else if(mode==1){
			benchQueueOnly(ff, threads);
		}else{
			benchFullParse(ff, threads);
		}
		t.stop("Total time: ");
	}

	private static void benchReadOnly(FileFormat ff){
		ByteFile bf=ByteFile.makeByteFile(ff);
		long lines=0, bytes=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			lines++;
			bytes+=line.length;
		}
		bf.close();
		System.err.println("Mode 0 (read only): "+lines+" lines, "+bytes+" bytes");
	}

	private static void benchQueueOnly(FileFormat ff, int threads){
		ByteFile bf=ByteFile.makeByteFile(ff);
		final JobQueue<ListNum<ArrayList<byte[]>>> queue=
				new JobQueue<>(threads+4, true, true, 0);
		final AtomicLong recordCount=new AtomicLong();

		Thread[] workers=new Thread[threads];
		for(int i=0; i<threads; i++){
			workers[i]=new Thread("Bench-"+i){
				@Override
				public void run(){
					long count=0;
					for(ListNum<ArrayList<byte[]>> batch=queue.take();
							batch!=null; batch=queue.take()){
						for(ArrayList<byte[]> record : batch){
							count++;
						}
					}
					recordCount.addAndGet(count);
				}
			};
			workers[i].start();
		}

		produce(bf, queue);

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}
		bf.close();
		System.err.println("Mode 1 (queue only): "+recordCount.get()+" records dispatched and discarded");
	}

	private static void benchFullParse(FileFormat ff, int threads){
		ByteFile bf=ByteFile.makeByteFile(ff);
		ConcurrentHashMap<Integer, Clade> map=new ConcurrentHashMap<>(16000);
		loadClades(bf, map, threads);
		bf.close();
		System.err.println("Mode 2 (full parse): "+map.size()+" clades loaded");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	static int loadThreads=-1;

	private static final int RECORDS_PER_BUNDLE=8;

}
