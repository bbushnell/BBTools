package jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ReadWrite;
import fun.DiskBench;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

public class TestFilesystem {
	
	public static void main(String[] args) throws InterruptedException{
		if(args.length>0){in=args[0];}
		if(args.length>1){out=args[1];}
		if(args.length>2){log=args[2];}
		if(args.length>3){size=Long.parseLong(args[3]);}
		if(args.length>4){ways=Integer.parseInt(args[4]);}
		if(args.length>5){interval_seconds=Long.parseLong(args[5]);}
		assert(args.length<=6);
		assert(interval_seconds>=0);
		assert(ways>0);

		ReadWrite.ZIPLEVEL=2;
		if("null".equalsIgnoreCase(out)){out=null;}
		if("null".equalsIgnoreCase(log)){log="stdout";}
		
		millis=interval_seconds*1000;

		Timer t=new Timer();
		for(int i=0; i<ways; i++){
			final String fname=in.replaceFirst("#", ""+i);
			File f=new File(fname);
			if(f.exists()){
				if(Math.abs(f.length()-size)>0.05*size){
					System.err.println("Existing file "+fname+" is the wrong size ("+f.length()+"); deleting.");
					f.delete();
				}
			}
			if(!f.exists()){
//				f.mkdirs();
				long bytes=DiskBench.writeRandomData(fname, size, t, true);
				System.err.println("Wrote "+fname+"; "+bytes+" bytes in "+t.timeInSeconds(3)+" seconds.");
			}
		}
		
		readQueue.put(new ByteBuilder(65536));
		readQueue.put(new ByteBuilder(65536));
		readQueue.put(new ByteBuilder(65536));
		
		runLoop();
	}

	private static void runLoop(){
		long nextTime=System.currentTimeMillis();
		final Timer t=new Timer();
		
		printLine(header());
		
		while(true){
			waitUntil(nextTime);
			final long time=System.currentTimeMillis();
			t.start();
			copy(in, out, iteration);
			t.stop();
			final long copyTime=t.elapsed;
			t.start();
			testMetadata();
			t.stop();
			final long metaTime=t.elapsed;
			t.start();
			testLs();
			t.stop();
			final long lsTime=t.elapsed;
			log(time, copyTime, metaTime, lsTime);
			nextTime+=millis;
			iteration=(iteration+1)%ways;
		}
	}
	
	private static void testMetadata(){
		{
			File f=new File("meta");
			if(!f.exists()){f.mkdir();}
		}
		if(fnames==null){
			fnames=new String[filesToCreate];
			for(int i=0; i<filesToCreate; i++){
				String fname="meta/"+i+".txt";
				fnames[i]=fname;
			}
		}
		for(String s : fnames){
			ReadWrite.write(s, s, false);
		}
		for(String s : fnames){
			InputStream is=ReadWrite.getRawInputStream(s, false);
			try {
				is.read(buffer);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		for(String s : fnames){
			new File(s).delete();
		}
	}
	
	private static int testLs(){
//		File dir=new File("x").getParentFile();
		File dir=new File(System.getProperty("user.dir"));
//		System.err.println(dir.getAbsolutePath());
		
		File[] array=dir.listFiles();
		assert(array!=null) : array;
		assert(array.length>0) : array.length;
		
		int count=0;
		for(File f : array){
			assert(f.canRead());
			count++;
		}
		return count;
	}
	
	private static void copy(String from0, final String to, int iter){
		final String from=from0.replaceFirst("#", ""+iter);
		{
			File f=new File(out);
			if(f.exists()){f.delete();}
		}
		final InputStream is=ReadWrite.getRawInputStream(from, false);
		WriteThread wt=new WriteThread(to);
		wt.start();
		for(int len=fillBuffer(is); len>0; len=fillBuffer(is)){
			ByteBuilder bb=readFetch();
			bb.append(buffer, 0, len);
			put(bb, wt.writeQueue);
		}
		put(poison, wt.writeQueue);
		try {
			is.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		while(wt.getState()!=Thread.State.TERMINATED){
			try {
				wt.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	private static int fillBuffer(final InputStream is){
		int len=0;
		while(len==0){
			try {
				len = is.read(buffer);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return len;
	}
	
	private static ByteBuilder readFetch(){
		ByteBuilder bb=null;
		while(bb==null){
			try {
				bb = readQueue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(bb.length()==0);
		assert(bb!=poison);
		return bb;
	}
	
	private static String header(){
		return "#time\tsize\tcopyTime\tMB/s\tmetaOps\tmetaTime\tops/s\tlsTime\tdate\n";
	}
	
	private static void log(final long time, final long copyTime, final long metaTime, final long lsTime){
		StringBuilder sb=new StringBuilder();
		sb.append(time).append('\t');
		sb.append(size).append('\t');
		sb.append(copyTime/1000000).append('\t');
		sb.append(Tools.format("%.2f", size*1000/Tools.max(1, (double)copyTime))).append('\t');
		sb.append(filesToCreate*3).append('\t');
		sb.append(metaTime/1000000).append('\t');
		sb.append(Tools.format("%.2f", (3*filesToCreate)*1000000000L/Tools.max(1, (double)metaTime))).append('\t');
		sb.append(Tools.format("%.2f", lsTime/1000000.0)).append('\t');
		sb.append(new Date(time)).append('\n');
		printLine(sb.toString());
	}
	
	private static void printLine(String line){
		if(log.equalsIgnoreCase("stdout") || log.equalsIgnoreCase("stdout.txt")){
			System.out.print(line);
		}else{
			ReadWrite.writeString(line, log, true);
		}
	}
	
	private static void waitUntil(final long time){
//		System.err.println("time="+time+", millis="+millis+", waiting until "+time);
		synchronized(waiter){
			long remaining=time-System.currentTimeMillis();
			while(remaining>9){//Minimum 10 millisecond wait.
				try {
					waiter.wait(remaining);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				remaining=time-System.currentTimeMillis();
			}
		}
	}
	
	private static void put(ByteBuilder bb, ArrayBlockingQueue<ByteBuilder> q) {
		while(bb!=null){
			try {
				q.put(bb);
				bb=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private static final class WriteThread extends Thread {
		
		WriteThread(String fname_){
			fname=fname_;
			outstream=ReadWrite.getRawOutputStream(fname, false, false);
		}
		
		@Override
		public void run(){
			ByteBuilder bb=writeFetch();
			while(bb!=poison){
				boolean success=write(bb);
				bb.setLength(0);
				put(bb, readQueue);
				bb=writeFetch();
			}
			try {
				outstream.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		private boolean write(ByteBuilder bb){
			assert(bb.length()>0);
			try {
				outstream.write(bb.array, 0, bb.length());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return false;
			}
			return true;
		}
		
		private ByteBuilder writeFetch(){
			ByteBuilder bb=null;
			while(bb==null){
				try {
					bb = writeQueue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return bb;
		}

		private final String fname;
		private final OutputStream outstream;
		
		public ArrayBlockingQueue<ByteBuilder> writeQueue=new ArrayBlockingQueue<ByteBuilder>(4);
		
	}
	
	
	private static int iteration;
	private static String in="foo#.txt";
	private static String out="bar.txt";
	private static String log="log.txt";
	private static int ways=40;
	private static int filesToCreate=1000;
	private static String[] fnames;
	private static long size=5000000000L;
	private static long interval_seconds=3600;
	private static long millis=interval_seconds*1000;
	private static final String waiter="waiter";
	private static final byte[] buffer=new byte[65536];
	private static final ByteBuilder poison=new ByteBuilder();
	private static ArrayBlockingQueue<ByteBuilder> readQueue=new ArrayBlockingQueue<ByteBuilder>(4);
	
}
