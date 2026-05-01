package clade;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import shared.Shared;
import shared.Tools;

/**
 * Selects threading strategy for clade loading based on file count
 * vs available threads. Ensures single-file inputs get full thread
 * utilization instead of falling back to single-threaded loading.
 *
 * @author Eru
 * @date April 30, 2026
 */
public class CladeLoaderAuto extends CladeObject {

	public ArrayList<Clade> loadFiles(ArrayList<String> in, boolean perContig,
			int minContig, long maxReads, boolean finish) {
		final int threads=Shared.threads();
		final int files=in.size();

		if(files==0){return new ArrayList<Clade>();}

		if(files==1){
			CladeLoaderSF loader=new CladeLoaderSF();
			ArrayList<Clade> result=loader.loadFile(in.get(0), perContig, minContig, maxReads, finish);
			readsProcessed+=loader.readsProcessed;
			basesProcessed+=loader.basesProcessed;
			return result;
		}

		if(files>=threads/2){
			CladeLoaderMF loader=new CladeLoaderMF();
			ArrayList<Clade> result=loader.loadFiles(in, perContig, minContig, maxReads, finish);
			readsProcessed+=loader.readsProcessed;
			basesProcessed+=loader.basesProcessed;
			return result;
		}

		return loadFewFiles(in, perContig, minContig, maxReads, finish, threads, files);
	}

	private ArrayList<Clade> loadFewFiles(ArrayList<String> in, boolean perContig,
			int minContig, long maxReads, boolean finish, int threads, int files) {
		final int threadsPerFile=Tools.max(1, threads/files);
		ExecutorService pool=Executors.newFixedThreadPool(files);
		@SuppressWarnings("unchecked")
		Future<ArrayList<Clade>>[] futures=new Future[files];
		final CladeLoaderSF[] loaders=new CladeLoaderSF[files];
		for(int i=0; i<files; i++){
			final String fname=in.get(i);
			final CladeLoaderSF loader=new CladeLoaderSF();
			loaders[i]=loader;
			futures[i]=pool.submit(()->
				loader.loadFile(fname, perContig, minContig, maxReads, finish, threadsPerFile)
			);
		}
		pool.shutdown();

		ArrayList<Clade> results=new ArrayList<Clade>();
		for(int i=0; i<files; i++){
			try{
				results.addAll(futures[i].get());
			}catch(Exception e){
				throw new RuntimeException(e);
			}
			readsProcessed+=loaders[i].readsProcessed;
			basesProcessed+=loaders[i].basesProcessed;
		}
		return results;
	}

	public long readsProcessed=0;
	public long basesProcessed=0;

}
