package clade;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import bin.AdjustEntropy;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * A/B benchmark + correctness check for Clade.parseClade (positional) vs Clade.parseCladeFlex
 * (order-independent, BytePrefixDispatcher).
 *
 * Correctness uses LITERAL canonical serialization: both Clades are re-serialized via toBytes() and the
 * bytes compared, so every parsed AND derived field is checked (gc/entropy/strandedness, full 16S/18S,
 * full DDL content -- not just presence).
 *
 * Two timings, both best-of-N and median (pass 0 = warmup, dropped):
 *   PARSE-ONLY: records pre-loaded in RAM, so it isolates parse+alloc.
 *   FULL COLD LOAD: read the .gz + decompress + group + parse + finish -- the operationally relevant number,
 *     since this parser runs once at DB-load startup. gzip/I-O is identical for both, so the delta is the
 *     parse cost within a realistic load.
 *
 * Run ON the cluster. Usage: CladeParseBench in=<spectra.gz> [passes=8]
 *
 * @author Noire
 */
public class CladeParseBench {

	public static void main(String[] args){
		String path=null;
		int passes=8;
		int threads=Shared.threads();
		for(String a : args){
			String[] s=a.split("=", 2);
			String k=s[0].toLowerCase(); String v=s.length>1 ? s[1] : null;
			if(k.equals("in") || k.equals("ref")){path=v;}
			else if(k.equals("passes")){passes=Integer.parseInt(v);}
			else if(k.equals("threads") || k.equals("t")){threads=Integer.parseInt(v);}
		}
		if(path==null){System.err.println("Usage: CladeParseBench in=<spectra.gz> [passes=8]"); return;}

		AdjustEntropy.load(4, 150);
		final ArrayList<ArrayList<byte[]>> records=loadRecords(path);
		final int n=records.size();
		System.err.println("Loaded "+n+" records into RAM.");

		final LineParser1 lp=new LineParser1('\t');

		//Correctness via canonical re-serialization.
		int mismatches=0, checked=0;
		for(ArrayList<byte[]> rec : records){
			if(rec.size()<=5){continue;}
			Clade a=Clade.parseClade(rec, lp);
			Clade b=Clade.parseCladeFlex(rec, lp);
			checked++;
			if(!serializedEqual(a, b)){
				mismatches++;
				if(mismatches<=5){System.err.println("MISMATCH tid="+a.taxID+" ("+a.name+")");}
			}
		}
		System.err.println("Correctness (canonical toBytes): "+(mismatches==0 ? "IDENTICAL" : mismatches+" MISMATCHES")+
			" over "+checked+" records.");
		if(mismatches>0){System.exit(1);}

		//PARSE-ONLY timing (in-RAM).
		final long[] oldP=new long[passes], flexP=new long[passes];
		for(int p=0; p<passes; p++){
			long sink=0;
			long t0=System.nanoTime();
			for(ArrayList<byte[]> rec : records){if(rec.size()>5){sink+=Clade.parseClade(rec, lp).taxID;}}
			long t1=System.nanoTime();
			for(ArrayList<byte[]> rec : records){if(rec.size()>5){sink+=Clade.parseCladeFlex(rec, lp).taxID;}}
			long t2=System.nanoTime();
			oldP[p]=t1-t0; flexP[p]=t2-t1;
			if(sink==0){System.err.println("(sink guard)");}
		}

		//FULL COLD LOAD timing (read .gz + parse + finish).
		final long[] oldL=new long[passes], flexL=new long[passes];
		for(int p=0; p<passes; p++){
			oldL[p]=fullLoad(path, false);
			flexL[p]=fullLoad(path, true);
		}

		report("PARSE-ONLY (in-RAM, single thread)", oldP, flexP, n);
		report("FULL COLD LOAD (.gz -> Clades, single thread)", oldL, flexL, n);

		//MULTITHREADED parse (in-RAM) -- the real load path is MT (CladeLoaderMT). Sweep thread counts.
		System.err.println("\n===== MULTITHREADED PARSE (in-RAM), best of 4 =====");
		System.err.println(String.format("%-8s %12s %12s %11s %10s", "threads", "old(ms)", "flex(ms)",
			"flex/old", "old Mr/s"));
		for(int T : threadSweep(threads)){
			long ob=Long.MAX_VALUE, fb=Long.MAX_VALUE;
			for(int p=0; p<4; p++){
				ob=Math.min(ob, parseMT(records, false, T));
				fb=Math.min(fb, parseMT(records, true, T));
			}
			System.err.println(String.format("%-8d %12.1f %12.1f %10.3fx %10.2f",
				T, ob/1e6, fb/1e6, (double)fb/ob, n/(ob/1e6)/1000));
		}
	}

	/** Parse every record across `threads` workers (work-stealing via an atomic index); return wall-clock ns. */
	static long parseMT(final ArrayList<ArrayList<byte[]>> records, final boolean flex, final int threads){
		final AtomicInteger idx=new AtomicInteger(0);
		final long[] sinks=new long[threads];
		final int n=records.size();
		final Thread[] ws=new Thread[threads];
		final long t0=System.nanoTime();
		for(int ti=0; ti<threads; ti++){
			final int tid=ti;
			ws[ti]=new Thread(){
				@Override public void run(){
					final LineParser1 lp=new LineParser1('\t');
					long sink=0; int i;
					while((i=idx.getAndIncrement())<n){
						ArrayList<byte[]> rec=records.get(i);
						if(rec.size()>5){sink+=(flex ? Clade.parseCladeFlex(rec, lp) : Clade.parseClade(rec, lp)).taxID;}
					}
					sinks[tid]=sink;
				}
			};
			ws[ti].start();
		}
		for(Thread w : ws){try{w.join();}catch(InterruptedException e){}}
		final long t=System.nanoTime()-t0;
		long total=0; for(long s : sinks){total+=s;}
		if(total==0){System.err.println("(mt sink guard)");}
		return t;
	}

	/** {1,2,4,8,16,...,threads} filtered to <=threads, distinct, ascending. */
	static int[] threadSweep(int threads){
		ArrayList<Integer> list=new ArrayList<Integer>();
		for(int t=1; t<threads; t*=2){list.add(t);}
		list.add(threads);
		int[] a=new int[list.size()];
		for(int i=0; i<a.length; i++){a[i]=list.get(i);}
		return a;
	}

	/** Read the .gz, group + parse + finish every record with the chosen parser; return elapsed ns. No cloning
	 *  (ByteFile returns fresh buffers, as CladeLoaderMT relies on) so this matches the real load path. */
	static long fullLoad(String path, boolean flex){
		final long t0=System.nanoTime();
		FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);
		ByteFile bf=ByteFile.makeByteFile(ff);
		final LineParser1 lp=new LineParser1('\t');
		ArrayList<byte[]> cur=new ArrayList<byte[]>(20);
		long sink=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(Tools.startsWith(line, '#') && cur.size()>5){
				sink+=(flex ? Clade.parseCladeFlex(cur, lp) : Clade.parseClade(cur, lp)).taxID;
				cur=new ArrayList<byte[]>(20);
				cur.add(line);
			}else{
				cur.add(line);
			}
		}
		if(cur.size()>5){sink+=(flex ? Clade.parseCladeFlex(cur, lp) : Clade.parseClade(cur, lp)).taxID;}
		bf.close();
		final long t=System.nanoTime()-t0;
		if(sink==0){System.err.println("(sink guard)");}
		return t;
	}

	static void report(String label, long[] oldT, long[] flexT, int n){
		//Drop pass 0 (warmup).
		final long ob=min(oldT, 1), fb=min(flexT, 1), om=median(oldT, 1), fm=median(flexT, 1);
		System.err.println("\n===== "+label+" ("+(oldT.length-1)+" timed passes) =====");
		System.err.println(String.format("parseClade     best=%.1fms median=%.1fms  (%.2f M rec/s best)",
			ob/1e6, om/1e6, n/(ob/1e6)/1000));
		System.err.println(String.format("parseCladeFlex best=%.1fms median=%.1fms  (%.2f M rec/s best)",
			fb/1e6, fm/1e6, n/(fb/1e6)/1000));
		System.err.println(String.format("flex/old       best=%.3fx  median=%.3fx  (%+.1f%% median)",
			(double)fb/ob, (double)fm/om, (fm-om)*100.0/om));
	}

	static long min(long[] a, int from){
		long m=Long.MAX_VALUE;
		for(int i=from; i<a.length; i++){m=Math.min(m, a[i]);}
		return m;
	}
	static long median(long[] a, int from){
		long[] c=Arrays.copyOfRange(a, from, a.length);
		Arrays.sort(c);
		return c[c.length/2];
	}

	/** Load records into RAM, grouping as CladeLoaderMT.produce does; owns the bytes (copyOf) so records persist. */
	static ArrayList<ArrayList<byte[]>> loadRecords(String path){
		FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);
		ByteFile bf=ByteFile.makeByteFile(ff);
		ArrayList<ArrayList<byte[]>> records=new ArrayList<ArrayList<byte[]>>(200000);
		ArrayList<byte[]> cur=new ArrayList<byte[]>(20);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			final byte[] copy=(line.length==0 ? line : Arrays.copyOf(line, line.length));
			if(Tools.startsWith(line, '#') && cur.size()>5){
				records.add(cur);
				cur=new ArrayList<byte[]>(20);
				cur.add(copy);
			}else{
				cur.add(copy);
			}
		}
		if(cur.size()>5){records.add(cur);}
		bf.close();
		return records;
	}

	/** Canonical-serialization equality: re-serialize both and compare bytes (covers every field + derived stat). */
	static boolean serializedEqual(Clade a, Clade b){
		ByteBuilder ba=a.toBytes(null), bb=b.toBytes(null);
		if(ba.length()!=bb.length()){return false;}
		return ba.toString().equals(bb.toString());
	}
}
