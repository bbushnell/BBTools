package idaligner;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ReadInputStream;
import structures.ByteBuilder;

public class Test {

	public static void main(String[] args) {

		args=new PreParser(args, null, false).args;
		long loops=(args.length>2 ? Integer.parseInt(args[2]) : 400);
		int threads=(args.length>3 ? Integer.parseInt(args[3]) : 1);
		int minLoops=(args.length>4 ? Integer.parseInt(args[4]) : 1);

		if (args.length>1) {
			String a=args[0], b=args[1];
			// Align provided sequences
			final byte[] seq1=toSequence(a);
			final byte[] seq2=toSequence(b);

			System.err.println(header());

			final long loops2=(seq2.length<500 ? loops : Tools.mid(1, loops, Tools.max(minLoops, threads)));
			test(new GlocalAligner(), seq1, seq2, loops2, threads);
			//	        test(new CrossCutAligner(), seq1, seq2, loops2, threads);
			test(new BandedAligner(), seq1, seq2, loops, threads);
			test(new DriftingAligner(), seq1, seq2, loops, threads);
			test(new WobbleAligner(), seq1, seq2, loops, threads);
			test(new QuantumAligner(), seq1, seq2, loops, threads);
			test(new QuabbleAligner(), seq1, seq2, loops, threads);
			test(new XDropHAligner(), seq1, seq2, loops, threads);
			test(new WaveFrontAligner2(), seq1, seq2, loops, threads);

			//	        test(new GlocalAligner(), seq1, seq2, loops2, threads);
			//	        test(new GlocalPlusAligner(), seq1, seq2, loops2, threads);
			////	        test(new GlocalPlusAligner2(), seq1, seq2, loops2, threads);
			////	        test(new GlocalPlusAligner3(), seq1, seq2, loops2, threads);
			////	        test(new GlocalPlusAligner4(), seq1, seq2, loops2, threads);
			//	        test(new GlocalPlusAligner5(), seq1, seq2, loops2, threads);
			//	        test(new CrossCutAligner(), seq1, seq2, loops2, threads);
			//	        test(new BandedAligner(), seq1, seq2, loops, threads);
			////	        test(new BandedAlignerM(), seq1, seq2, loops, threads);
			//	        test(new BandedPlusAligner(), seq1, seq2, loops, threads);
			//	        test(new BandedPlusAligner2(), seq1, seq2, loops, threads);
			////	        test(new BandedPlusAlignerInt(), seq1, seq2, loops, threads);
			//	        test(new BandedByteAligner(), seq1, seq2, loops, threads);
			//	        test(new DriftingAligner(), seq1, seq2, loops, threads);
			////	        test(new DriftingAlignerM(), seq1, seq2, loops, threads);
			//	        test(new DriftingPlusAligner(), seq1, seq2, loops, threads);
			//	        test(new DriftingPlusAligner2(), seq1, seq2, loops, threads);
			//	        test(new DriftingPlusAligner3(), seq1, seq2, loops, threads);//Broken?
			//	        test(new WobbleAligner(), seq1, seq2, loops, threads);
			////	        test(new WobblePlusAligner(), seq1, seq2, loops, threads);
			////	        test(new WobblePlusAligner2(), seq1, seq2, loops, threads);
			//	        test(new WobblePlusAligner3(), seq1, seq2, loops, threads);
			//	        test(new WobblePlusAligner5(), seq1, seq2, loops, threads);//Broken?
			//	        test(new QuantumAligner(), seq1, seq2, loops, threads);
			////	        test(new QuantumPlusAligner(), seq1, seq2, loops, threads);
			////	        test(new QuantumAlignerM(), seq1, seq2, loops, threads);
			////	        test(new QuantumPlusAligner3(), seq1, seq2, loops, threads);
			//	        test(new QuantumPlusAligner4(), seq1, seq2, loops, threads);
			//	        test(new SquabbleAligner(), seq1, seq2, loops, threads);
			////	        test(new WaveFrontAligner(), seq1, seq2, loops, threads);
			////	        test(new SingleStateAlignerFlat2(), seq1, seq2, loops, threads);
			//	        

		}
	}

	public static boolean validate(IDAligner ida) {
		float f;

		f=test(ida, "A", "A");
		assert(f==1f);

		f=test(ida, "T", "A");
		assert(f==0.0f);

		f=test(ida, "AA", "AA");
		assert(f==1f);

		f=test(ida, "AAA", "A");
		assert(f==0.33333333f);

		f=test(ida, "CCC", "A");
		assert(f==0f);

		f=test(ida, "AA", "AGA");
		assert(f==0.6666667f) : f;

		f=test(ida, "AGA", "AA");
		assert(f==0.6666667f) : f;

		f=test(ida, "AT", "AA");
		assert(f==0.5f) : f;

		f=test(ida, "AAA", "AAA");
		assert(f==1f);

		f=test(ida, "AAAT", "AAAA");
		assert(f==0.75f);

		f=test(ida, "ACGA", "AAAA");
		assert(f==0.5f);

		f=test(ida, "AAAA", "AAAAA");
		assert(f==1f);

		f=test(ida, "AAGAA", "AAAA");
		assert(f==0.8f);

		f=test(ida, "AAAA", "AAGAA");
		assert(f==0.8f) : f;

		f=test(ida, "AAAAAA", "AAAAAA");
		assert(f==1f) : f;

		f=test(ida, "CCCCCC", "AAAAAA");
		assert(f==0f) : f;

		f=test(ida, "AAAAAA", "AAAAAA");
		assert(f==1f) : f;

		f=test(ida, "AAAAAAA", "AAAAAAA");
		assert(f==1f) : f;

		f=test(ida, "AAATAAA", "AAAAAAA");
		assert(f==6f/7) : f;

		f=test(ida, "AAAAAAAA", "AAAAAAAA");
		assert(f==1f) : f+"!="+1;

		f=test(ida, "AAATAAAA", "AAAAAAAA");
		assert(f==7f/8) : f;

		f=test(ida, "AAATAAAAAA", "AAAAAAAAAA");
		assert(f==9f/10) : f;

		f=test(ida, "AAAAAAAAAAAA", "AAAAAAAAAAAA");
		assert(f==1f) : f;

		f=test(ida, "AAATAAAAAAAA", "AAAAAAAAAAAA");
		assert(f==11f/12) : f;

		f=test(ida, "AAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAA");
		assert(f==1f) : f;

		f=test(ida, "CCCCCCCCCCCCCCCC", "A");
		assert(f==0f) : f;

		f=test(ida, "AAAAAATTTTAAAAAA", "AAAAAAAAAAAAAAAA");
		assert(f==0.75f) : f;

		f=test(ida, "AAAAAAAAAAAA", "AAAAAATTTTAAAAAA");
		assert(f==0.75f) : f;

		return true;
	}

	public static final void print(long[] curr, String name) {
		System.err.print(name+" Score:\t");
		for(int i=0; i<curr.length; i++) {System.err.print((curr[i]>>42)+" ");}
		System.err.print("\n"+name+" Dels: \t");
		for(int i=0; i<curr.length; i++) {System.err.print((((curr[i]>>21)&0xFFFF)+" "));}
		System.err.println();
	}

	public static byte[] toSequence(String a) {
		if(a.length()<100 && new File(a).isFile()) {
			return ReadInputStream.toReads(a, FileFormat.FA, 1).get(0).bases;
		}else {return a.getBytes();}
	}

	public static String header() {
		return "Name     \tANI\trStart\trStop\tLoops\tSpace%\tTime";
	}

	public static void printResults(IDAligner ida, byte[] a, byte[] b, float id, int[] pos, 
			long iters, int threads, Timer t) {
		AlignmentStats stats=new AlignmentStats();
		stats.setAndSolve(pos, a.length, b.length);
		printResults(ida, a, b, stats, iters, threads, t);
	}

	public static void printResults(IDAligner ida, byte[] a, byte[] b, 
			AlignmentStats stats, long iters, int threads, Timer t) {
		long loops=ida.loops()/iters;
		t.stop();
		String time=t.timeInSeconds(3);

		float stateSpace=a.length*b.length;
		float fraction=loops/stateSpace;
		ByteBuilder bb=new ByteBuilder();

		bb.append(ida.name());
		while(bb.length()<9) {bb.append(' ');}
		bb.tab();
		bb.appendt(stats.identity, 4);
		bb.appendt(stats.rStart);
		bb.appendt(stats.rStop);
		bb.appendt(loops);
		bb.appendt(fraction*100, 3);
		bb.append(time);
		System.err.println(bb);
	}

	public static void testAndPrint(String name, String[] args) {
		IDAligner ida=createNewInstance(name);
		testAndPrint(ida, args);
	}

	public static <C extends IDAligner> void testAndPrint(Class<C> c, String[] args) {
		IDAligner ida=createNewInstance(c);
		testAndPrint(ida, args);
	}

	public static <C extends IDAligner> void testAndPrint(IDAligner ida, String[] args) {
		PrintStream outstream=System.err;
		String output=null;
		String query=null, ref=null;
		boolean debug=false;
		boolean verbose=false;
		boolean printOps=false;
		boolean global=false;
		boolean validate=false;
		boolean statsmode=false;
		String aligner=null;
		int iters=1;

//		Parser parser=new Parser();
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			//Custom args
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				//TODO
			}else if(a.equals("stats") || a.equals("alignmentstats")){
				statsmode=Parse.parseBoolean(b);
			}else if(a.equals("pos") || a.equals("posvector")){
				statsmode=!Parse.parseBoolean(b);
			}else if(a.equals("trace") || a.equals("traceback")){
				DO_TRACE=Parse.parseBoolean(b);
			}else if(a.equals("printtrace")){
				PRINT_TRACE=Parse.parseBoolean(b);
			}else if(a.equals("query") || a.equals("q") || a.equals("in") || a.equals("in1")){
				query=b;
			}else if(a.equals("ref") || a.equals("r") || a.equals("in2")){
				ref=b;
			}else if(a.equals("iters") || a.equals("iterations")){
				iters=Integer.parseInt(b);
			}else if(a.equals("printops")){
				printOps=Tracer.PRINT_OPS=Parse.parseBoolean(b);
			}else if(a.equals("debug")){
				debug=Parse.parseBoolean(b);
			}else if(a.equals("global")){
				global=Tracer.GLOBAL=Parse.parseBoolean(b);
			}else if(a.equals("test") || a.equals("validate")){
				validate=Parse.parseBoolean(b);
			}else if(a.equals("type") || a.equals("aligner")){
				aligner=arg;
			}

			//Parser args
			else if(Parser.parseStatic(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}

			//Positional args
			else if(i==0 && (Tools.looksLikeInputSequenceStream(arg) || arg.indexOf('=')<0)){
				query=arg;
			}else if(i==1 && (Tools.looksLikeInputSequenceStream(arg) || arg.indexOf('=')<0)){
				ref=arg;
			}else if(i==2 && ("null".equalsIgnoreCase(arg) || Tools.looksLikeOutputStream(arg))) {
				output=("null".equalsIgnoreCase(arg) ? null : arg);
			}else if(i==3 && Tools.isNumeric(arg)) {
				iters=Integer.parseInt(arg);
			}

			//Unknown args
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		DO_TRACE|=PRINT_TRACE;
		assert((ida==null) != (aligner==null));
		if(ida==null) {ida=createNewInstance(aligner);}
		if(validate || (query==null & ref==null)) {
			Test.validate(ida);
			return;
		}

		final Class<? extends IDAligner> c=ida.getClass();
		if(debug) {setOrCatch(c, "debug", debug);}
		if(global) {setOrCatch(c, "GLOBAL", global);}
		if(printOps) {setOrCatch(c, "PRINT_OPS", global);}
		if(output!=null) {setOrCatch(c, "output", output);}
		if(statsmode){testAndPrintStats(ida, args[0], args[1], iters);}
		else{testAndPrintPos(ida, args[0], args[1], iters);}
	}

	private static boolean setOrCatch(Class<?> c, String name, Object value) {
		try {
			Field f = c.getDeclaredField(name);
			if(java.lang.reflect.Modifier.isStatic(f.getModifiers())) {
				//	            f.setAccessible(true);// Should be public
				f.set(null, value); // null for static field
			}
		} catch(NoSuchFieldException | IllegalAccessException e) {
			System.err.println("Can't set field "+name+" for class "+c.getSimpleName());
			return false;
		}
		return true;
	}

	public static void testAndPrintPos(IDAligner ida, String a, String b, int iters) {
		int[] pos=new int[4];
		Timer t=new Timer();
		float id=0;
		ida.setLoops(0);
		byte[] a2=toSequence(a);
		byte[] b2=toSequence(b);
		ByteBuilder bb=new ByteBuilder();
		bb.append(ida.name()).colon();
		while(bb.length()<9) {bb.append(' ');}
		for(int i=0; i<iters; i++) {id=ida.align(a2, b2, pos);}
		long loops=ida.loops();
		float percent=loops*100f/(iters*1f*a2.length*b2.length);
		bb.append("\tid=").append(id, 5).append("\tcoords=(").append(pos[0]).comma();
		bb.append(pos[1]).append(")\tloops=").append(loops/iters).tab().append(percent, 2).percent();
		t.stopAndStart(bb.tab().toString());
	}

	public static void testAndPrintStats(IDAligner ida, String a, String b, int iters) {
		AlignmentStats stats=new AlignmentStats();
		stats.doTrace=DO_TRACE;
		Timer t=new Timer();
		float id=0;
		ida.setLoops(0);
		byte[] a2=toSequence(a);
		byte[] b2=toSequence(b);
		ByteBuilder bb=new ByteBuilder();
		bb.append(ida.name()).colon();
		while(bb.length()<9) {bb.append(' ');}
		for(int i=0; i<iters; i++) {id=ida.align(a2, b2, stats);}
		long loops=ida.loops();
		float percent=loops*100f/(iters*1f*a2.length*b2.length);
		bb.append("\tid=").append(id, 5).append("\tcoords=(").append(stats.rStart).comma();
		bb.append(stats.rStop).append(")\tloops=").append(loops/iters).tab().append(percent, 2).percent();
		t.stopAndStart(bb.tab().toString());
		assert(!PRINT_TRACE || stats.matchString!=null) : 
			DO_TRACE+", "+PRINT_TRACE+", "+stats.doTrace+", "+stats.matchString;
		if(PRINT_TRACE && stats.matchString!=null) {
			t.outstream.println("\n"+new String(stats.matchString));
		}
	}

	public static float test(IDAligner ida, String a, String b) {
		return testST(ida, toSequence(a), toSequence(b), 1, 0);
	}

	public static float test(IDAligner ida, String a, String b, int iters) {
		return testST(ida, toSequence(a), toSequence(b), iters, 0);
	}

	public static float test(IDAligner ida, byte[] a, byte[] b, long maxIters, int threads) {
		if(threads<2) {return testST(ida, a, b, maxIters, 1000);}
		else {return testMT(ida, a, b, maxIters, threads);}
	}

	public static float testST(IDAligner ida, byte[] a, byte[] b, long maxIters, long sleepTime) {
		if(sleepTime>0) {
			try {Thread.sleep(sleepTime);}//Cool the CPU 
			catch (InterruptedException e) {e.printStackTrace();}
		}
		ida.setLoops(0);

		Timer t=new Timer();
		float id=0;
		int[] pos=new int[2];
		for(int i=0; i<maxIters; i++) {id=ida.align(a, b, pos);}
		t.stop();

		printResults(ida, a, b, id, pos, maxIters, 1, t);
		return id;
	}

	@SuppressWarnings("unchecked")
	public static <T> T createNewInstance(T existingObject) {
		try {
			Class<?> c=existingObject.getClass();
			return (T)c.getDeclaredConstructor().newInstance();
		}catch(ReflectiveOperationException e){
			throw new RuntimeException(e);
		}
	}

	@SuppressWarnings("unchecked")
	public static <T> T createNewInstance(Class<T> c) {
		try {
			return (T)c.getDeclaredConstructor().newInstance();
		}catch(ReflectiveOperationException e){
			throw new RuntimeException(e);
		}
	}

	@SuppressWarnings("unchecked")
	public static IDAligner createNewInstance(String name) {
		Class<? extends IDAligner> c=getClass(name);
		return createNewInstance(c);
	}
	
	@SuppressWarnings("unchecked")
	private static Class<? extends IDAligner> getClass(String name){
		if(!name.startsWith("aligner.")) {name="aligner."+name;}
		try{return (Class<? extends IDAligner>)Class.forName(name);}
		catch(Exception e) {throw new RuntimeException(e);}
	}

	public static float testMT(IDAligner ida, byte[] a, byte[] b, long maxIters, int threads) {
		try {Thread.sleep(1000);}//Cool the CPU 
		catch (InterruptedException e) {e.printStackTrace();}

		if(threads<1) {
			threads=Math.max(1, (int)Math.min(Math.min(maxIters, 256), Shared.threads()));
		}

		//Disable loop tracking in multithreaded version
		if(threads>1) {ida.setLoops(-1);}

		Timer t=new Timer();
		AtomicLong iters=new AtomicLong(0);
		ArrayList<Runner> list=new ArrayList<Runner>(threads);
		for(int i=0; i<threads; i++) {
			Runner runner=new Runner(a, b, createNewInstance(ida), iters, maxIters, i);
			list.add(runner);
		}

		//    	ExecutorService pool=Executors.newFixedThreadPool(threads);
		//    	for(Runner runner : list){
		//    	    pool.submit(runner);
		//    	}
		//    	pool.shutdown();
		//    	boolean completed=false;
		//    	while(!completed){
		//    	    try{
		//    	        completed=pool.awaitTermination(60, TimeUnit.SECONDS);
		//    	    }catch(InterruptedException e){
		//    	        Thread.currentThread().interrupt();
		//    	        e.printStackTrace();
		//    	    }
		//    	}

		boolean success=processThreads(list, threads);

		float id=-1;
		int[] pos=null;
		for(Runner x : list) {
			synchronized(x) {
				if(id<x.id) {
					id=x.id;
					pos=x.pos;
				}
			}
		}

		printResults(ida, a, b, id, pos, maxIters, threads, t);
		return id;
	}

	/** This uses ThreadPools as an experiment, but they are very slow. */
	public static boolean processThreads(List<? extends Runnable> list, int threads) {
		// Create pool once
		ExecutorService pool=Executors.newFixedThreadPool(threads);

		boolean success=processThreadsInPool(list, pool);

		// Only call shutdown when completely done with the pool
		pool.shutdown();

		return success;
	}

	// Multiple uses of the same pool
	public static boolean processThreadsInPool(List<? extends Runnable> list, ExecutorService pool) {
		// Submit tasks
		List<Future<?>> futures=new ArrayList<>();
		for(Runnable runner : list){
			futures.add(pool.submit(runner));
		}

		// Wait for all tasks to complete
		// Check for exceptions
		boolean success=true;
		for(Future<?> future : futures){
			try{
				future.get();
			}catch(Exception e){
				success=false;
				e.printStackTrace();
			}
		}

		if(!success){
			System.err.println("Warning: Some alignment threads failed");
		}
		return success;
	}

	private static class Runner implements Runnable{

		Runner(byte[] q, byte[] r, IDAligner ida_, AtomicLong iters_, long maxIters_, int tnum_){
			query=q;
			ref=r;
			ida=ida_;
			iters=iters_;
			maxIters=maxIters_;
			tnum=tnum_;
		}

		@Override
		public void run() {
			synchronized(this) {
				//Make a local copy
				pos=new int[2];
				query=Arrays.copyOf(query, query.length);
				ref=Arrays.copyOf(ref, ref.length);

				//Loop
				for(long iter=iters.getAndIncrement(); iter<maxIters; iter=iters.getAndIncrement()) {
					id=ida.align(query, ref, pos);
				}
			}
		}

		byte[] query;
		byte[] ref;
		final IDAligner ida;
		final AtomicLong iters;
		final long maxIters;
		final int tnum;
		int[] pos;
		float id=-1;

	}
	
	static boolean DO_TRACE=false;
	static boolean PRINT_TRACE=false;

}
