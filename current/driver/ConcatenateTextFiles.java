package driver;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import tracker.ReadStats;

/**
 * Utility for concatenating multiple text files into a single output file.
 * Supports both individual files and directory processing with buffered I/O.
 * Uses multi-threaded writing to handle large file operations efficiently.
 * @author Brian Bushnell
 */
public class ConcatenateTextFiles {
	
	/** Format: infile1,infile2,...infileN,outfile */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		if(ReadWrite.ZIPLEVEL<6){ReadWrite.ZIPLEVEL=6;}
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else{
				concatenate(args[i].split(","));
			}
				
		}
		t.stop();
		System.out.println();
		System.out.println("Time: \t"+t);
	}

	/**
	 * Concatenates input files into a single output file using buffered processing.
	 * Creates a write thread for output and processes each input file sequentially.
	 * @param split Array of file paths where last element is the output file
	 */
	private static void concatenate(String[] split) {
		String outname=split[split.length-1];
		assert(overwrite || !new File(outname).exists()) : outname+" exists.";
		
		WriteThread wt=new WriteThread(outname);
		wt.start();
		
		
		ArrayList<String>[] bufferptr=new ArrayList[] {new ArrayList<String>(LIST_SIZE)};
		
		for(int i=0; i<split.length-1; i++){
			processTerm(split[i], bufferptr, wt);
		}
		
		ArrayList<String> buffer=bufferptr[0];
		if(buffer==null){
			wt.add(new ArrayList<String>(1));
		}else if(buffer.isEmpty()){
			wt.add(buffer);
		}else{
			wt.add(buffer);
			wt.add(new ArrayList<String>(1));
		}
			
	}
	
	/**
	 * Processes a single input term (file or directory) for concatenation.
	 * If term is a file, reads lines into buffer. If directory, recursively processes contents.
	 * Sends full buffers to the write thread to maintain memory efficiency.
	 *
	 * @param term File path or directory path to process
	 * @param bufferptr Array containing current buffer reference
	 * @param wt WriteThread instance for output handling
	 */
	private static void processTerm(String term, ArrayList<String>[] bufferptr, WriteThread wt){
		
		System.out.println("Processing term "+term);
		
		File f=new File(term);
		if(!f.isDirectory()){

			TextFile tf=new TextFile(term, false);
			
			ArrayList<String> buffer=bufferptr[0];
			
			String s=null;
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				buffer.add(s);

				//			System.out.println("Added to buffer");
				if(buffer.size()>=LIST_SIZE){
					//				System.out.println("Sent buffer");
					
//					System.out.println("****** "+term+" ******");
//					for(String b : buffer){
//						System.out.println(b);
//					}
					
					wt.add(buffer);
					bufferptr[0]=buffer=new ArrayList<String>(LIST_SIZE);
				}
			}
			tf.close();
		}else{
			assert(f.isDirectory());
			File[] contents=f.listFiles();
			for(File c : contents){
				String abs=c.getAbsolutePath();
				if(!abs.equals(wt.fname)){
//					System.out.println(c+" == "+new File(wt.fname)+" : "+c.equals(new File(wt.fname)));
					processTerm(abs, bufferptr, wt);
				}
			}
		}
	}
	
	/**
	 * Dedicated thread for writing concatenated text data to output file.
	 * Uses blocking queue to receive string buffers from processing threads.
	 * Handles file I/O operations asynchronously for better performance.
	 */
	private static class WriteThread extends Thread{
		
		/**
		 * Creates a WriteThread for the specified output file.
		 * Initializes output stream and print writer for file operations.
		 * @param fname_ Output file path
		 */
		public WriteThread(String fname_){
			String temp=fname_;
			try {
				temp=new File(fname_).getCanonicalPath();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			fname=temp;
			os=ReadWrite.getOutputStream(fname, append, true, true);
			writer=new PrintWriter(os);
		}
		
		/**
		 * Adds a string buffer to the write queue for output processing.
		 * Blocks until queue space is available to prevent memory overflow.
		 * @param list Buffer containing text lines to write
		 */
		public void add(ArrayList<String> list){
			assert(list!=null);
			while(list!=null){
//				System.out.println("Adding list to queue "+queue.size());
				try {
					queue.put(list);
					list=null;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		@Override
		public void run(){
			
			ArrayList<String> list=null;
			while(list==null){
//				System.out.println("Waiting for list...");
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
//				System.out.println("Took list of size "+(list==null ? "null" : list.size()+""));
				if(list!=null){
					if(list.isEmpty()){
						ReadWrite.finishWriting(writer, os, fname, allowSubprocess);
						return;
					}
					for(String s : list){
						if(s!=null){writer.println(s);}
					}
				}
				list=null;
			}
		}
		
		/** Output stream for writing to the target file */
		private final OutputStream os;
		/** Print writer for formatted text output */
		private final PrintWriter writer;
		/**
		 * Blocking queue for thread-safe buffer transfer between processing and writing
		 */
		private final ArrayBlockingQueue<ArrayList<String>> queue=new ArrayBlockingQueue<ArrayList<String>>(MAX_LISTS);
		/** Canonical path of the output file */
		final String fname;
		
	}

	/** Maximum number of string buffers allowed in the write queue */
	public static final int MAX_LISTS=8;
	/** Number of text lines stored in each buffer before sending to write thread */
	public static final int LIST_SIZE=100;
	/** Whether to overwrite existing output files */
	public static boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	public static boolean append=false;
	/** Whether to allow subprocess execution for file operations */
	public static boolean allowSubprocess=true;
	
}
